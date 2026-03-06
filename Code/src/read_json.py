"""Read and preprocess AML mutation trees from JSON format.

This script:
1. Loads mutation trees from a JSON file (aml_mutation_trees.js).
2. Merges duplicate child nodes within each tree.
3. Removes redundant ancestor-descendant label duplicates.
4. Assigns numeric IDs to mutations (using Jaccard similarity > 0.5
   to match multi-label nodes to existing entries).
5. Exports the processed tree map as a pickle file and mutation/sample
   mappings as CSV files.
"""

import json
import os
import re
import pickle
from typing import Dict, List, Optional, Set

import pandas as pd
from anytree import PreOrderIter
from anytree.importer import JsonImporter

# Jaccard similarity threshold for matching multi-label nodes to existing mutations
JACCARD_THRESHOLD = 0.5


def check_duplicate_labels(root) -> List:
    """Find duplicate mutation labels in a tree.

    Args:
        root: The root node of an anytree tree.

    Returns:
        List of label values that appear more than once in the tree.
    """
    labels = [node.label for node in PreOrderIter(root)]
    duplicates = set(label for label in labels if labels.count(label) > 1)
    return list(duplicates)


def read_json_trees(path: str) -> Dict:
    """Load mutation tree data from a JSON file.

    Args:
        path: Path to the JSON file.

    Returns:
        Dictionary mapping sample names to their tree data.
    """
    with open(path) as json_file:
        return json.load(json_file)


def extract_aml_number(sample_name: str) -> Optional[str]:
    """Extract the base AML patient ID from a sample name.

    Example: 'AML-42-001' -> 'AML-42'

    Args:
        sample_name: Full sample identifier string.

    Returns:
        Base patient ID, or None if the pattern doesn't match.
    """
    match = re.match(r'(AML-\d+)-\d+', sample_name)
    if match:
        return match.group(1)
    return None


def merge_duplicate_children(root) -> None:
    """Merge children with duplicate labels within a tree.

    When two children of the same parent share a label, the second child's
    subtree is grafted onto the first and the duplicate is removed.
    Also removes nodes whose label matches an ancestor's label.

    Args:
        root: The root node of the tree (modified in place).
    """
    for node in PreOrderIter(root):
        # Merge children with the same label
        seen_labels: Dict = {}
        for child in node.children:
            if child.label not in seen_labels:
                seen_labels[child.label] = child
            else:
                if child.children != ():
                    existing_node = seen_labels[child.label]
                    existing_node.children = existing_node.children + child.children
                child.parent = None

        # Remove nodes whose label duplicates an ancestor's label
        for ancestor in node.ancestors:
            if ancestor.label == node.label:
                if ancestor.children != ():
                    parent_node = node.parent
                    for child in node.children:
                        parent_node.children = parent_node.children + (child,)
                node.parent = None


def aggregate_multilabel_nodes(tree_map: Dict, parallel_trees: Dict, mutation_index: Dict) -> int:
    """Aggregate ancestor labels into multi-label nodes for trees with duplicates.

    For nodes whose label appears in the duplicate list, the node's label
    is extended to include ancestor labels, creating a multi-label node.

    Args:
        tree_map: Dictionary mapping patient IDs to their tree dicts.
        parallel_trees: Dict of {patient: {tree_id: [duplicate_labels]}}.
        mutation_index: Running dict mapping label strings to numeric IDs.

    Returns:
        Next available mutation index.
    """
    next_id = max(mutation_index.values()) + 1 if mutation_index else 0

    for patient_id in tree_map:
        for tree_id in tree_map[patient_id]:
            tree_root = tree_map[patient_id][tree_id]
            for node in PreOrderIter(tree_root):
                # Expand labels for nodes in trees with duplicates
                if (patient_id in parallel_trees
                        and tree_id in parallel_trees[patient_id]
                        and node.label in parallel_trees[patient_id][tree_id]):
                    label = node.label
                    node.label = [label]
                    for ancestor in node.ancestors:
                        if ancestor.label != 'Root':
                            if isinstance(ancestor.label, list):
                                node.label += ancestor.label
                            else:
                                node.label += [ancestor.label]
                    # Collapse single-element lists back to scalar
                    if len(node.label) == 1:
                        node.label = node.label[0]

                # Assign numeric ID to mutation label
                label_key = str(node.label)
                if label_key not in mutation_index:
                    if isinstance(node.label, list):
                        # Try to match via Jaccard similarity to existing labels
                        node_label_set = set(node.label)
                        primary_label = node.label[0]
                        matched = False
                        for existing_key in mutation_index:
                            existing_set = set()
                            if existing_key.startswith('['):
                                existing_labels = existing_key.strip('[]').split(',')
                                existing_primary = existing_labels[0]
                                existing_set = set(existing_labels)
                            else:
                                existing_primary = existing_key
                                existing_set.add(existing_key)
                            jaccard = len(node_label_set & existing_set) / len(node_label_set | existing_set)
                            if jaccard > JACCARD_THRESHOLD and primary_label == existing_primary:
                                mutation_index[label_key] = mutation_index[existing_key]
                                matched = True
                                break
                        if not matched:
                            mutation_index[label_key] = next_id
                            next_id += 1
                    else:
                        mutation_index[label_key] = next_id
                        next_id += 1

                node.name = mutation_index[label_key]

                if node.label == 'Root':
                    node.data = str(tree_id)
                else:
                    node.data = node.node_id

    return next_id


def process_aml_trees(input_path: str, output_dir: str) -> None:
    """Full pipeline: load JSON trees, preprocess, and export results.

    Args:
        input_path: Path to the input JSON file (aml_mutation_trees.js).
        output_dir: Directory for output files (tree_map.pickle, CSVs).
    """
    sample_map = read_json_trees(input_path)
    tree_map: Dict = {}
    importer = JsonImporter()
    tree_counter = 200  # Starting ID for tree instances
    mutation_index: Dict[str, int] = {}
    sample_to_number: Dict[str, List[int]] = {}

    # --- Phase 1: Import and merge duplicates ---
    for sample_name, sample_data in sample_map.items():
        root = importer.import_(json.dumps(sample_data["event_tree"]))
        merge_duplicate_children(root)

        sample_to_number[sample_name] = [tree_counter]
        base_id = extract_aml_number(sample_name)
        patient_key = base_id + '-001'
        if patient_key not in tree_map:
            tree_map[patient_key] = {}
        tree_map[patient_key][tree_counter] = root
        tree_counter += 1

    # --- Phase 2: Detect duplicate labels across trees ---
    parallel_trees: Dict = {}
    for patient_id in tree_map:
        for tree_id in tree_map[patient_id]:
            duplicate_labels = check_duplicate_labels(tree_map[patient_id][tree_id])
            if len(duplicate_labels) > 0:
                if patient_id not in parallel_trees:
                    parallel_trees[patient_id] = {}
                parallel_trees[patient_id][tree_id] = duplicate_labels

    # --- Phase 3: Aggregate multi-label nodes and assign mutation IDs ---
    aggregate_multilabel_nodes(tree_map, parallel_trees, mutation_index)

    # --- Phase 4: Export results ---
    mutation_df = pd.DataFrame.from_dict(mutation_index, orient='index', columns=['Index'])
    mutation_df.to_csv(os.path.join(output_dir, 'aml_mut_map.csv'))

    sample_df = pd.DataFrame.from_dict(sample_to_number, orient='index', columns=['Index'])
    sample_df.to_csv(os.path.join(output_dir, 'aml_sample_map.csv'))

    sorted_tree_map = dict(sorted(tree_map.items(), key=lambda kv: int(kv[0].split('-')[1])))

    with open(os.path.join(output_dir, 'tree_map.pickle'), 'wb') as handle:
        pickle.dump(sorted_tree_map, handle, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    current_directory = os.getcwd()
    process_aml_trees(
        input_path=os.path.join(current_directory, "AML", "aml_mutation_trees.js"),
        output_dir=os.path.join(current_directory, "AML"),
    )
