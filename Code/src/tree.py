"""Tree data structure for representing cancer mutation trees.

Wraps AnyTree's NodeMixin to provide a typed tree node that stores
mutation labels and arbitrary payload data.
"""

from typing import Generic, Iterable, Optional, TypeVar

import anytree

NameType = TypeVar("NameType", int, str)
DataType = TypeVar("DataType")


class TreeNode(Generic[NameType, DataType], anytree.NodeMixin):
    """Tree node compatible with AnyTree.

    Can be annotated with the label type and the data type.

    Attributes:
        name: Identifier of the node (int or str).
        data: Arbitrary payload data (e.g., node index or metadata).
        mutation_label: The mutation label(s) assigned to this node.
        parent: The parent node, or None for the root.
        children: Iterable of child nodes.
    """

    def __init__(
        self,
        name: NameType,
        data: DataType = None,
        mutation_label: NameType = None,
        parent: Optional["TreeNode"] = None,
        children: Optional[Iterable["TreeNode"]] = None,
    ) -> None:
        """
        Args:
            name: Identifier of the node (compatible with declared type).
            data: Payload data attached to this node.
            mutation_label: The mutation label for this node (e.g., gene name).
            parent: The parent node in the tree.
            children: Iterable of child nodes.
        """
        self.name = name
        self.data = data
        self.mutation_label = mutation_label
        self.parent = parent
        if children:
            self.children = children

    def __str__(self) -> str:
        """Render the tree as an indented string."""
        ret = ""
        for pre, _, node in anytree.RenderTree(self):
            ret += f"{pre}{node.name}\n"
        return ret


__all__ = ["TreeNode", "DataType"]
