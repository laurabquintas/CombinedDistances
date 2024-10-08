"""This submodule defines our data format for trees, based
on AnyTree.
"""
from typing import Generic, Iterable, Optional, TypeVar

import anytree

NameType = TypeVar("NameType", int, str)
DataType = TypeVar("DataType")


class TreeNode(Generic[NameType, DataType], anytree.NodeMixin):
    """Tree node compatible with AnyTree.
    It can be annotated with the label type and the data type.
    Attrs:
        name: identifier of the node (compatible with declared type)
        data: payload data
        parent: the node of the parent
        children: iterable with the children nodes
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
            name: identifier of the node (compatible with declared type)
            data: payload data
            parent: the node of the parent
            children: iterable with the children nodes
        """
        self.name = name
        self.data = data
        self.mutation_label = mutation_label
        self.parent = parent
        if children:
            self.children = children

    def __str__(self) -> str:
        """Casts to str."""
        ret = ""
        for pre, _, node in anytree.RenderTree(self):
            ret += f"{pre}{node.name}\n"
        return ret


__all__ = ["TreeNode", "DataType"]