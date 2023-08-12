#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import networkx as nx
import argparse
from collections import defaultdict

ITER_LIMIT_PER_LOCALMOVE = -1
MIN = 0.0000001

def generate_network(community_size, num_communities, intra_edges, inter_edges, p1, p2, seed=None):
    if seed is not None:
        np.random.seed(seed)

    G_positive = nx.Graph()
    G_negative = nx.Graph()

    for i in range(num_communities):
        G_tmp = nx.gnm_random_graph(community_size, intra_edges)
        G_tmp = nx.relabel_nodes(G_tmp, {node: node + i*community_size for node in G_tmp.nodes()})
        G_positive = nx.compose(G_positive, G_tmp)

    inter_edge_candidates = [(i, j) for i in range(community_size*num_communities) for j in range(i+1, community_size*num_communities) if abs(i//community_size - j//community_size) == 1]
    inter_edge_selected = np.random.choice(len(inter_edge_candidates), inter_edges * num_communities, replace=False)
    negative_edges = [inter_edge_candidates[i] for i in inter_edge_selected]
    G_negative.add_edges_from(negative_edges)

    G_positive.add_nodes_from(G_negative.nodes())
    G_negative.add_nodes_from(G_positive.nodes())

    positive_edges = list(G_positive.edges())
    for edge_index in np.random.choice(len(positive_edges), int(intra_edges * p1), replace=False):
        edge = positive_edges[edge_index]
        if edge[0]//community_size == edge[1]//community_size and G_positive.degree(edge[0]) > 1 and G_positive.degree(edge[1]) > 1:
            G_positive.remove_edge(*edge)
            G_negative.add_edge(*edge)

    negative_edges = list(G_negative.edges())
    for edge_index in np.random.choice(len(negative_edges), int(inter_edges * p2), replace=False):
        edge = negative_edges[edge_index]
        if edge[0]//community_size != edge[1]//community_size and G_negative.degree(edge[0]) > 1 and G_negative.degree(edge[1]) > 1:
            G_negative.remove_edge(*edge)
            G_positive.add_edge(*edge)

    return G_positive, G_negative

def renumber(dictionary):
    values = set(dictionary.values())
    renumbering = dict(zip(values, range(len(values))))
    return {k: renumbering[v] for k, v in dictionary.items()}

class GraphInfo:
    def __init__(self, ref_graph, graph, partition, weight_key):
        self.graph_whole = nx.Graph()
        self.graph_whole.add_nodes_from(ref_graph.nodes(data=True))
#        self.graph_whole.add_edges_from(graph.edges(data=True))
        for u, v, data in graph.edges(data=True):
            if u in self.graph_whole.nodes() and v in self.graph_whole.nodes():
                self.graph_whole.add_edge(u, v, **data)

        self.weight_key = weight_key
        self.node_degrees = dict()
        self.loops = dict()
        self.size = float(self.graph_whole.size(weight=weight_key))
        self.community_degrees = dict()
        self.community_sizes = dict()

        for node in self.graph_whole.nodes():
            node_degree = float(self.graph_whole.degree(node, weight=weight_key))
            self.node_degrees[node] = node_degree
            loop_edge = graph.get_edge_data(node, node, default={weight_key: 0.})
            self.loops[node] = float(loop_edge.get(weight_key))

        self._initialize_community_info(partition)

    def _initialize_community_info(self, partition):
        self.community_degrees = {community: 0 for community in set(partition.values())}
        self.community_sizes   = {community: 0 for community in set(partition.values())}

        for node, community in partition.items():
            node_degree = self.node_degrees[node]
            loop = self.loops[node]
            self.community_degrees[community] += node_degree
            self.community_sizes[community] += loop

    def calc_modularity(self, partition):
        m = self.size
        modularity = 0.
        for community in set(partition.values()):
            m_c = self.community_sizes.get(community, 0.)
            K_c = self.community_degrees.get(community, 0.)
            modularity += m_c / m -  ((K_c / (2. * m)) ** 2)
        return modularity

    def _remove_node(self, node, partition, linksum_dict):
        community = partition.get(node)
        linksum = linksum_dict.get(community, 0.)
        self.community_degrees[community] = self.community_degrees.get(community, 0.) - self.node_degrees.get(node, 0.)
        self.community_sizes[community] = self.community_sizes.get(community, 0.) - linksum - self.loops.get(node, 0.)

    def _insert_node(self, node, partition, linksum_dict, community):
        linksum = linksum_dict.get(community, 0.)
        self.community_degrees[community] = self.community_degrees.get(community, 0.) + self.node_degrees.get(node, 0.)
        self.community_sizes[community] = self.community_sizes.get(community, 0.) + linksum + self.loops.get(node, 0.)

    def get_linksum_dict(self, node, partition):
        weight_key = self.weight_key
        graph = self.graph_whole
        linksum_dict = defaultdict(float)

        for neighbor_node, edge in graph[node].items():
            if neighbor_node != node:
                neighbor_community = partition[neighbor_node]
                linksum_dict[neighbor_community] += edge.get(weight_key, 0)

        return linksum_dict

    def _delta_q_1(self, node, linksum_dict, partition):
        community = partition.get(node)
        ki = self.node_degrees.get(node, 0.)
        ac2m = self.community_degrees.get(community, 0.)
        m = self.size
        linksum = linksum_dict.get(community, 0.)  # linksum_dictにはそのnodeのエッジ重みの和が格納されている
        result = - linksum + ac2m*ki/(2*m) - ki**2/(2*m)
        return result

    def _delta_q_2(self, node, linksum_dict, neighboring_community):
        ki = self.node_degrees.get(node, 0.)
        ac2m = self.community_degrees.get(neighboring_community, 0.)
        m = self.size
        linksum = linksum_dict.get(neighboring_community, 0.)
        result = linksum - ac2m*ki/(2*m)
        return result

class LouvainSigned:
    def __init__(self, positive_graph, negative_graph, alpha, random_generator = np.random.default_rng(), seed=None, mode="positive"):
        if not isinstance(alpha, (int, float)):
            print("Error: Alpha value must be a number. Please provide a valid value.")
            return
        elif not 0 <= alpha <= 1:
            print("Error: Alpha value must be in the range [0, 1]. Please provide a valid value.")
            return

        self.weight_key = 'weight'
        if seed != None:
            self.random_generator = np.random.default_rng(seed=seed)
        else:
            self.random_generator = random_generator
        self.graph_whole = self.get_signed_graph(positive_graph, negative_graph, mode)

        self.partition = {node: i for i, node in enumerate(self.graph_whole.nodes())}
#        print(self.partition)

        self.ginfo_pos = GraphInfo(self.graph_whole, positive_graph, self.partition, self.weight_key)
        self.ginfo_neg = GraphInfo(self.graph_whole, negative_graph, self.partition, self.weight_key)

        self.dendrogram = list()
        self.alpha = alpha
        self.mode = mode

    def get_signed_graph(self, positive_graph, negative_graph, mode):
        graph = nx.Graph()
        if mode == "Full":
            graph.add_nodes_from(positive_graph.nodes(data=True))
            graph.add_nodes_from(negative_graph.nodes(data=True))
            # graph.add_edges_from(positive_graph.edges(data=True), sign='positive')
            # graph.add_edges_from(negative_graph.edges(data=True), sign='negative')
        else:
            graph.add_nodes_from(positive_graph.nodes(data=True))

        return graph

    def _randomize(self, items):
        randomized_items = list(items)
        self.random_generator.shuffle(randomized_items)
        return randomized_items

    def _modularity(self, partition):
        alpha = self.alpha
        Q_pos = self.ginfo_pos.calc_modularity(partition)
        Q_neg = self.ginfo_neg.calc_modularity(partition)
#        print ("Q_pos, Q_neg", Q_pos, Q_neg, alpha * Q_pos + (1-alpha) * (1 - Q_neg))
        return alpha * Q_pos + (1-alpha) * (1 - Q_neg)

    def _move_nodes(self):
        modified = True
        nb_pass_done = 0
        new_Q = self._modularity(self.partition)
        alpha = self.alpha

        while modified and nb_pass_done != ITER_LIMIT_PER_LOCALMOVE:
            current_Q = new_Q
            modified = False
            nb_pass_done += 1

            for node in self._randomize(self.graph_whole.nodes()):
                original_community = self.partition.get(node)
                linksum_dict_pos = self.ginfo_pos.get_linksum_dict(node, self.partition)
                linksum_dict_neg = self.ginfo_neg.get_linksum_dict(node, self.partition)
                q1_pos = self.ginfo_pos._delta_q_1(node, linksum_dict_pos, self.partition)
                q1_neg = self.ginfo_neg._delta_q_1(node, linksum_dict_neg, self.partition)
                self.ginfo_pos._remove_node(node, self.partition, linksum_dict_pos)
                self.ginfo_neg._remove_node(node, self.partition, linksum_dict_neg)
                self.partition[node] = -1

                best_community = original_community
                best_increase = 0
                for neighboring_community in self._randomize(linksum_dict_pos.keys()):
                    q2_pos = self.ginfo_pos._delta_q_2(node, linksum_dict_pos, neighboring_community)
                    q2_neg = self.ginfo_neg._delta_q_2(node, linksum_dict_neg, neighboring_community)

                    delta_Q = alpha * (q1_pos + q2_pos) - (1 - alpha) * (1 - q1_neg - q2_neg)
#                    print ("deltaQ", neighboring_community, delta_Q)
                    if delta_Q > best_increase:
                        best_increase = delta_Q
                        best_community = neighboring_community

                self.ginfo_pos._insert_node(node, self.partition, linksum_dict_pos, best_community)
                self.ginfo_neg._insert_node(node, self.partition, linksum_dict_neg, best_community)
                self.partition[node] = best_community
#                print ("best com", best_community)
                if best_community != original_community:
                    modified = True

            new_Q = self._modularity(self.partition)
            if new_Q - current_Q < MIN:
                break

    def _aggregate_nodes(self, partition, graph):
        weight = self.weight_key
        aggregated_graph = nx.Graph()
        aggregated_graph.add_nodes_from(partition.values())

        for node1, node2, _data in graph.edges(data=True):
            edge_weight = _data.get(weight, 1)
            c1 = partition[node1]
            c2 = partition[node2]
            w_prec = aggregated_graph.get_edge_data(c1, c2, {weight: 0}).get(weight, 1)
            aggregated_graph.add_edge(c1, c2, **{weight: w_prec + edge_weight})
        return aggregated_graph

    def generate_dendrogram(self, mode):
        current_graph_pos = self.ginfo_pos.graph_whole.copy()
        current_graph_neg = self.ginfo_neg.graph_whole.copy()
        partition_list = list()
        Q = -1.0

        while True:
            self._move_nodes()
            renumberd_partition = renumber(self.partition)
            partition_list.append(renumberd_partition)
            current_graph_pos = self._aggregate_nodes(renumberd_partition, current_graph_pos)
            current_graph_neg = self._aggregate_nodes(renumberd_partition, current_graph_neg)
            self.graph_whole = self.get_signed_graph(current_graph_pos, current_graph_neg, mode)

            self.partition = {node: i for i, node in enumerate(self.graph_whole.nodes())}
            self.ginfo_pos = GraphInfo(self.graph_whole, current_graph_pos, self.partition, self.weight_key)
            self.ginfo_neg = GraphInfo(self.graph_whole, current_graph_neg, self.partition, self.weight_key)
            new_Q = self._modularity(self.partition)
#            print("generate_dendrogram", new_Q)
            if new_Q - Q < MIN:
                break
            Q = new_Q

        self.dendrogram = partition_list[:]

    def best_partition(self):
        self.generate_dendrogram(self.mode)
        partition = self.dendrogram[0].copy()
        for level in range(1, len(self.dendrogram)):
            for node, community in partition.items():
                partition[node] = self.dendrogram[level][community]
        return partition

def main():
    parser = argparse.ArgumentParser(prog='eeisp')
    parser.add_argument("matrix", help="Input matrix", type=str)
    parser.add_argument("output", help="Output prefix", type=str)
    parser.add_argument("--alpha", help="alpha parameter (from 0 to 1, default: 0.5)", type=float, default=0.5)
    parser.add_argument("--seed", help="seed for randamization", type=int)
    parser.add_argument("-v", "--version", action='version', version='%(prog)s 0.6.0')

    args = parser.parse_args()
    print(args)

    alpha = args.alpha
    community_size = 8
    num_communities = 4
    intra_edges = 12
    inter_edges = 10
    p1 = 0.05
    p2 = 0.05

    G_positive, G_negative = generate_network(community_size, num_communities, intra_edges, inter_edges, p1, p2)

    print(f"G_positive: Nodes={G_positive.number_of_nodes()}, Edges={G_positive.number_of_edges()}")
    print(f"G_negative: Nodes={G_negative.number_of_nodes()}, Edges={G_negative.number_of_edges()}")
    print(f"G: Nodes={nx.compose(G_positive, G_negative).number_of_nodes()}, Edges={nx.compose(G_positive, G_negative).number_of_edges()}")

    l = LouvainSigned(G_positive, G_negative, alpha, seed=10)
    partition = l.best_partition()
    print(f"partition: {partition}")

if __name__ == "__main__":
    main()
