import sys
import os
sys.path.insert(0, os.path.abspath('.'))
"""
de Bruijn Graph Genome Assembler
Implements k-mer graph construction, Eulerian path finding, and contig assembly.
"""

from collections import defaultdict, deque
from typing import List, Dict, Optional
import networkx as nx


class DeBruijnAssembler:
    """
    de Bruijn graph assembler for short-read genome assembly.

    Algorithm:
        1. Break reads into k-mers
        2. Build directed graph: (k-1)-mer nodes, k-mer edges
        3. Find Eulerian path through graph
        4. Reconstruct sequence from path
    """

    def __init__(self, k: int = 21):
        self.k = k
        self.graph = defaultdict(list)   # adjacency list
        self.in_degree = defaultdict(int)
        self.out_degree = defaultdict(int)

    def add_read(self, read: str):
        """Decompose a read into k-mers and add edges to graph."""
        read = read.upper().strip()
        if len(read) < self.k:
            return
        for i in range(len(read) - self.k + 1):
            kmer = read[i:i+self.k]
            prefix = kmer[:-1]   # (k-1)-mer left node
            suffix = kmer[1:]    # (k-1)-mer right node
            self.graph[prefix].append(suffix)
            self.out_degree[prefix] += 1
            self.in_degree[suffix]  += 1
            # Ensure all nodes appear
            if prefix not in self.in_degree:  self.in_degree[prefix]  = 0
            if suffix not in self.out_degree: self.out_degree[suffix] = 0

    def add_reads(self, reads: List[str]):
        """Add multiple reads to the graph."""
        for read in reads:
            self.add_read(read)

    def find_start_node(self) -> Optional[str]:
        """
        Find Eulerian path start: node with out_degree - in_degree = 1.
        If none, return any node with outgoing edges.
        """
        for node in self.graph:
            if self.out_degree[node] - self.in_degree[node] == 1:
                return node
        # Eulerian circuit: start anywhere
        return next(iter(self.graph), None)

    def find_eulerian_path(self) -> List[str]:
        """Find Eulerian path using Hierholzer's algorithm."""
        graph_copy = {k: list(v) for k, v in self.graph.items()}
        start = self.find_start_node()
        if start is None:
            return []

        stack  = [start]
        path   = []

        while stack:
            node = stack[-1]
            if graph_copy.get(node):
                next_node = graph_copy[node].pop(0)
                stack.append(next_node)
            else:
                path.append(stack.pop())

        return list(reversed(path))

    def path_to_sequence(self, path: List[str]) -> str:
        """Reconstruct sequence from Eulerian path."""
        if not path:
            return ''
        seq = path[0]
        for node in path[1:]:
            seq += node[-1]
        return seq

    def assemble(self, reads: List[str]) -> str:
        """Full pipeline: reads → sequence."""
        self.add_reads(reads)
        path = self.find_eulerian_path()
        return self.path_to_sequence(path)

    def get_contigs(self, reads: List[str]) -> List[str]:
        """
        Generate contigs by finding all non-branching paths.
        More realistic than single Eulerian path for real data.
        """
        self.add_reads(reads)
        visited_edges = set()
        contigs = []

        for start in list(self.graph.keys()):
            # Start a contig at branching points or unvisited nodes
            if self.in_degree[start] != 1 or self.out_degree[start] != 1:
                for neighbor in self.graph[start]:
                    edge = (start, neighbor)
                    if edge not in visited_edges:
                        contig_nodes = [start, neighbor]
                        visited_edges.add(edge)
                        # Extend while non-branching
                        while (self.in_degree[neighbor] == 1 and
                               self.out_degree[neighbor] == 1):
                            next_n = self.graph[neighbor][0]
                            e = (neighbor, next_n)
                            if e in visited_edges: break
                            visited_edges.add(e)
                            contig_nodes.append(next_n)
                            neighbor = next_n
                        contigs.append(self.path_to_sequence(contig_nodes))
        return contigs
