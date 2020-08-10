import pickle
import matplotlib.pyplot as plt
import pygraphviz as pvg
import networkx as nx
import pandas as pd
import numpy as np
import requests
import mygene
import json
import time
import sys
import os
from Bio import SeqIO
from copy import deepcopy
from py2cytoscape import cyrest
from biothings_client import get_client
from networkx.drawing.nx_agraph import write_dot
from compute_impact import ImpactFactor
from networkx.algorithms import community
from networkx.algorithms import components

from bokeh.io import show
from bokeh.models import BasicTicker, ColorBar, LinearColorMapper, PrintfTickFormatter
from bokeh.plotting import figure
from bokeh.plotting import figure, output_file, save


def ensemble2symbol_dict(ids):
    mg = mygene.MyGeneInfo()
    unique_ids = np.unique(ids)
    responses = mg.getgenes(unique_ids, fields='symbol')
    print(responses[:3])
    return {x: x['symbol'] for x in responses if 'symbol' in x}


def load_interactome():
    res_fname = 'data/interactome.tsv'
    data = None
    if not os.path.isfile(res_fname):
        fname = 'data/HI-II-14.tsv'
        data = pd.read_csv(fname, sep='\t', header=None).values
        symbols_dict = ensemble2symbol_dict(data.flatten())
        # print(symbols)
        n, m = np.shape(data)
        for i in range(n):
            for j in range(m):
                if data[i, j] in symbols_dict:
                    data[i, j] = symbols_dict[data[i, j]]
                else:
                    data[i, j] = None
        data = pd.DataFrame(data)
        data.dropna(inplace=True)
        data.to_csv(res_fname, header=None, index=None, sep='\t')
    else:
        data = pd.read_csv(res_fname, sep='\t', header=None)
    return data


def load_disease_modules():
    # fname = 'data/all_gene_disease_associations.tsv'
    fname = 'data/curated_gene_disease_associations.tsv'
    data = pd.read_csv(fname, sep='\t')[['geneSymbol', 'diseaseName']].values
    modules = {}
    for i in range(np.shape(data)[0]):
        disease = data[i][1]
        if disease in modules:
            modules[disease].append(data[i][0])
        else:
            modules[disease] = [data[i][0]]
    return modules


def preprocess_if_data(data):
    data = data[['isoform_name', 'bs conservation', 'nucleotides per bs']]
    data.columns = ['isoform', 'bs', 'bs_size']

    def parse_bs(x):
        return np.array(x.strip('[] ').split(', '), dtype=float)

    def parse_bs_size(x):
        x = x.strip('[] ').split(', ')
        return np.array([y.strip('[]') for y in x], dtype=int)

    data['bs'] = data['bs'].apply(parse_bs)
    data['bs_size'] = data['bs_size'].apply(parse_bs_size)
    return data


def load_gtex():
    out_fname = 'data/gtex.pickle'
    fname = 'data/tissues_SMTS_tpm_mean_curated.tsv'
    if os.path.isfile(out_fname):
        return pickle.load(open(out_fname, 'rb'))
    data = pd.read_csv(fname, sep=',', index_col=0)
    trans_list_fname = 'data/transcripts_list.pickle'
    trans_dict_fname = 'data/transcripts_dict.pickle'
    transcripts_dict = None
    if os.path.isfile(trans_dict_fname):
        transcripts_dict = pickle.load(open(trans_dict_fname, 'rb'))
        transcripts_list = pickle.load(open(trans_list_fname, 'rb'))
    else:
        transcripts_list, transcripts_dict, no_result = convert_genes(data.index, conv_from='ensembl.transcript', conv_to='symbol', species=9606)
        pickle.dump(transcripts_list, open(trans_list_fname, 'wb'))
        pickle.dump(transcripts_dict, open(trans_dict_fname, 'wb'))
        # print(no_result)
        # print(len(no_result))
    gene_list = transcripts_list  # [transcripts_dict[x].upper() if x in transcripts_dict else None for x in data.index]
    data['symbol'] = gene_list
    data.fillna(value=0, inplace=True)
    pickle.dump(data, open(out_fname, 'wb'))

    return data


def load_tissue_level():
    out_fname = 'data/tissue_expression.pickle'
    if os.path.isfile(out_fname):
        return pickle.load(open(out_fname, 'rb'))
    gtex = load_gtex()
    pickle.dump(gtex, open(out_fname, 'wb'))
    return gtex


def load_transcript_level():
    out_fname = 'data/transcript_expression.pickle'
    if os.path.isfile(out_fname):
        return pickle.load(open(out_fname, 'rb'))
    gtex = load_gtex()
    gene_symbol = gtex['symbol']
    expression = gtex.sum(axis=1, numeric_only=True)
    transcript_expression = {'expression': expression, 'symbol': gene_symbol}
    transcript_expression = pd.DataFrame(transcript_expression)
    pickle.dump(transcript_expression, open(out_fname, 'wb'))
    return transcript_expression


def load_gene_level():
    out_fname = 'data/gene_expression.pickle'
    if os.path.isfile(out_fname):
        return pickle.load(open(out_fname, 'rb'))
    transcript_expression = load_transcript_level()
    gene_expression = transcript_expression.groupby(['symbol']).sum()
    pickle.dump(gene_expression, open(out_fname, 'wb'))
    return gene_expression


def compute_tissue_if():
    gtex = load_tissue_level()
    gtex = gtex.loc[~gtex.index.duplicated(keep='first')]
    # print(gtex.index, len(gtex.index), len(np.unique(gtex.index)))
    gtex_dict = gtex.to_dict(orient='index')
    return gtex
    # transcripts =


# comm = """
def compute_gene_if(expressions, bs, imp):
    # expression = expression.drop([0])
    isoforms = bs['isoform']
    groups = set([x.split('-')[0] for x in isoforms])
    impact_factor = {}
    for group in groups:
        if group not in expressions.index:
            continue
        idx = [i for i, s in enumerate(isoforms) if group == s.split('-')[0]]
        bs_change = np.concatenate(bs['bs'][idx].values)
        impact_factor[group] = imp.gene_if(expressions.loc[group]['expression'], -1, -1, bs_change)
    return impact_factor
# """


def compute_mock_gene_if(expression, bs, imp):
    isoforms = bs['isoform']
    groups = set([x.split('-')[0] for x in isoforms])
    impact_factor = {}
    for group in groups:
        if group not in expression.index:
            continue
        idx = [i for i, s in enumerate(isoforms) if group == s.split('-')[0]]
        bs_change = np.concatenate(bs['bs'][idx].values)
        impact_factor[group] = imp.gene_if(1, -1, -1, bs_change)
    return impact_factor


def load_impact_factor():
    out_fname = 'data/IF.pickle'
    fname = 'data/BS_Summary_combined.tsv'
    tissue_if_fname = 'data/tissue_IF.pickle'
    if os.path.isfile(out_fname):
        return pickle.load(open(out_fname, 'rb'))
    impact_factor = {}
    gene_expression = load_gene_level()
    data = pd.read_csv(fname, sep='\t')
    data = preprocess_if_data(data)
    imp = ImpactFactor(strategy='average')

    mapping = map_compas_ids()
    mapping = {v: k for k, v in mapping.items()}

    transcript_expression = load_transcript_level()

    tissue_if = None
    if os.path.isfile(tissue_if_fname):
        tissue_if = pickle.load(open(tissue_if_fname, 'rb'))
    else:
        tissue_expression = load_tissue_level()

        tissue_expression.index = [mapping[x] if x in mapping else None for x in tissue_expression.index]
        tissue_expression = tissue_expression.loc[tissue_expression.index.dropna()]
        data.index = [x.strip() for x in data['isoform']]

        for iso in tissue_expression.index:
            for col in tissue_expression.columns[:-1]:
                if iso in data.index:
                    expr = tissue_expression.loc[iso, col]
                    if np.shape(np.array([expr]).flatten())[0] > 1:
                        expr = expr.values[0]
                    bs = np.array([data.loc[iso, 'bs']]).flatten()
                    imp_factor = imp.gene_if(expr, -1, -1, bs)
                    tissue_expression.loc[iso, col] = imp_factor

        tissue_if = tissue_expression
        pickle.dump(tissue_expression, open(tissue_if_fname, 'wb'))

    genes_dataset = deepcopy(tissue_if)
    # print(tissue_if.mean(axis=0))
    # print(tissue_if.std(axis=0))

    genes_if_proto = deepcopy(genes_dataset)
    genes_if_proto.index = genes_dataset['symbol']
    genes_if_proto.drop(['symbol'], axis=1, inplace=True)

    genes_if_tissue_max = deepcopy(genes_if_proto.groupby(['symbol']).max())
    genes_if_tissue_max = genes_if_tissue_max.loc[genes_if_tissue_max.index.drop(0)]
    impact_factor['genes_tissue_max'] = deepcopy(genes_if_tissue_max)

    genes_if_tissue_var = deepcopy(genes_if_proto.groupby(['symbol']).std())
    genes_if_tissue_var = genes_if_tissue_var.loc[genes_if_tissue_var.index.drop(0)]
    impact_factor['genes_tissue_var'] = genes_if_tissue_var

    genes_if_max = genes_if_tissue_max.max(axis=1)
    genes_if_max = genes_if_max.to_dict()

    genes_var = deepcopy(genes_dataset)
    genes_var.index = genes_dataset['symbol']
    genes_var.drop(['symbol'], axis=1, inplace=True)
    genes_var.drop([0], axis=0, inplace=True)
    genes_if_var = {}
    for sym in genes_var.index:
        ifs = genes_var.loc[sym].values
        genes_if_var[sym] = np.std(ifs.flatten())

    impact_factor['genes_max'] = genes_if_max
    impact_factor['genes_var'] = genes_if_var
    pickle.dump(impact_factor, open(out_fname, 'wb'))

    return impact_factor


def remove_duplicate_positive_edges(data):
    data = data.sort_values('predictions')
    return data.drop_duplicates(data.columns[[2, 4]])


def color_nodes(G, nodes, color):
    cutoff = 0.9
    for node in G.node:
        if node in nodes:
            G.node[node]['color'] = color
        if node in nodes and color == 'red':
            if G.node[node]['impact_factor'] > cutoff:
                G.node[node]['color'] = 'red'
            else:
                G.node[node]['color'] = 'yellow'
        else:
            if G.node[node]['impact_factor'] > cutoff:
                G.node[node]['color'] = 'turquoise'
            else:
                G.node[node]['color'] = 'grey'
        # print(node, G.node[node]['impact_factor'], G.node[node]['color'])
    return G


def is_hub(G, node_label, hub_threshold):
    return len(G[node_label]) >= hub_threshold


def prune_network(G, node_attrs, node_attr_values, edge_attrs, edge_attr_values, hub_threshold=5, affect_neighbors=True):
    nodes2remain = []
    for node1, node2, edge in G.edges.data():
        for attr, values in zip(edge_attrs, edge_attr_values):
            if edge[attr] in values:
                nodes2remain += [node1, node2]
                break
        for attr, values in zip(node_attrs, node_attr_values):
            if affect_neighbors:
                if (G.node[node1][attr] in values or G.node[node2][attr] in values) or (is_hub(G, node1, hub_threshold) or is_hub(G, node2, hub_threshold)):
                    nodes2remain += [node1, node2]
                    break
            else:
                if is_hub(G, node1, hub_threshold) or is_hub(G, node2, hub_threshold):
                    nodes2remain += [node1, node2]
                elif G.node[node1][attr] in values:
                    nodes2remain.append(node1)
                elif G.node[node2][attr] in values:
                    nodes2remain.append(node2)

    nodes2remove = np.array(G.nodes)[[x not in nodes2remain for x in G.nodes]]
    G.remove_nodes_from(nodes2remove)
    return G


# PyGraphViz


def default_node_style(A):
    A.node_attr['style'] = 'filled'
    A.node_attr['shape'] = 'circle'
    A.node_attr['fixed_size'] = 'true'
    return A

# Cytoscape


def create_new_network(json_data):
    PORT_NUMBER = 1234
    BASE = 'http://localhost:' + str(PORT_NUMBER) + '/v1/'

    # Header for posting data to the server as JSON
    HEADERS = {'Content-Type': 'application/json'}

    res = requests.post(BASE + 'networks?collection=ppi',
                        data=json.dumps(json_data), headers=HEADERS)
    network_id = res.json()['networkSUID']
    return network_id, res


def update_cytoscape_style(style_name, mappings):
    cytoscape.vizmap.update_style(title=style_name, mappings=mappings)
    cytoscape.vizmap.apply(styles=style_name)


def apply_diabetes_style(cytoscape, statistics):
    style_name = 'T2D'
    defaults_dic = {"NODE_SHAPE": "circle",
                    "NODE_SIZE": "30",
                    "NODE_FILL_COLOR": "#024b4d",
                    "NODE_LABEL_POSITION": 'S,NW,c,0.00,3.00',
                    "NODE_LABEL_FONT_SIZE": "14",
                    "EDGE_TRANSPARENCY": "255"}
    defaults_list = cytoscape.vizmap.simple_defaults(defaults_dic)
    cytoscape.vizmap.create_style(title=style_name, defaults=defaults_list)
    time.sleep(2)
    mappings = []

    node_color_map = cytoscape.vizmap.mapVisualProperty(
        visualProperty='NODE_FILL_COLOR', mappingType='discrete',
        mappingColumn='color', discrete=[['#420D09', '#800000', '#FA8072', '#1890F0', '#A8D8F0', '#ebf4f5'], ['#420D09', '#800000', '#FA8072', '#1890F0', '#A8D8F0', '#ebf4f5']])
    # mappingColumn='color', discrete=[['turquoise', 'red', 'yellow', 'grey'], ['#024b4d', '#911d23', '#f7ca18', '#d8d7de']])
    # mappingColumn='color', discrete=[['grey', 'red'], ['#d8d7de', '#911d23']])
    mappings.append(node_color_map)

    node_label_map = cytoscape.vizmap.mapVisualProperty(visualProperty='NODE_LABEL', mappingType='passthrough', mappingColumn='name')
    mappings.append(node_label_map)

    node_shape_map = cytoscape.vizmap.mapVisualProperty(
        visualProperty='NODE_SHAPE', mappingType='discrete', mappingColumn='shape', discrete=[['circle', 'triangle', 'hexagon'], ['circle', 'diamond', 'diamond']])
    update_cytoscape_style(style_name, [node_shape_map])
    mappings.append(node_shape_map)

    edge_color_map = cytoscape.vizmap.mapVisualProperty(visualProperty='EDGE_STROKE_UNSELECTED_PAINT', table='edge', mappingType='discrete', mappingColumn='predictions', discrete=[
        ['-100', '-1', '0', '1'], ['#000000', '#cf001e', '#7ccdf2', '#7ccdf2']])    # grey '#a5a7a8' blue 7ccdf2 cyan cceeff
    mappings.append(edge_color_map)

    edge_visibility_map = cytoscape.vizmap.mapVisualProperty(visualProperty='EDGE_VISIBLE', table='edge', mappingType='discrete', mappingColumn='predictions', discrete=[
        ['-100', '-1', '0', '1'], ['false', 'true', 'true', 'true']])    # grey '#a5a7a8' blue 7ccdf2 cyan cceeff
    mappings.append(edge_visibility_map)

    # edge_width_map = cytoscape.vizmap.mapVisualProperty(
    #    visualProperty='EDGE_WIDTH', table='edge', mappingType='discrete', mappingColumn='predictions', discrete=[['-1', '0', '1'], ['5', '0.7', '1.5']])
    # edge_width_map = cytoscape.vizmap.mapVisualProperty(
    #    visualProperty='EDGE_WIDTH', table='edge', mappingType='continuous', mappingColumn='betweenness', discrete=[['-1', '0', '1'], ['5', '0.7', '1.5']])
    # mappings.append(edge_width_map)

    if 'betweenness' in statistics:
        # edge_transparency_map_continuous = cytoscape.vizmap.mapVisualProperty(
        #    visualProperty='EDGE_TRANSPARENCY', table='edge', mappingType='continuous', mappingColumn='betweenness', lower=[betweenness[0], 80], center=[betweenness[1], 180], upper=[betweenness[2], 255])
        # print('Continuous mapping!')
        # mappings.append(edge_transparency_map_continuous)
        betweenness = statistics['betweenness'].values()
        betweenness = [min(betweenness), 0.5 *
                       (min(betweenness)+max(betweenness)), max(betweenness)]

        edge_width_map_continuous = cytoscape.vizmap.mapVisualProperty(
            visualProperty='EDGE_WIDTH', table='edge', mappingType='continuous', mappingColumn='betweenness', lower=[betweenness[0], 0.5], center=[betweenness[1], 2], upper=[betweenness[2], 10])
        mappings.append(edge_width_map_continuous)

    if 'transparency' in statistics:
        transparency = statistics['transparency']
        # str_transp = [str(x) for x in np.unique(transparency)]
        # edge_transparency_map_discreet = cytoscape.vizmap.mapVisualProperty(
        #    visualProperty='EDGE_TRANSPARENCY', table='edge', mappingType='discrete', mappingColumn='transparency', discrete=[str_transp, str_transp])
        # mappings.append(edge_transparency_map_discreet)

    if 'border_width' in statistics:
        border_width = statistics['border_width']
        str_width = [str(x) for x in np.unique(border_width)]
        border_width_map = cytoscape.vizmap.mapVisualProperty(
            visualProperty='NODE_BORDER_WIDTH', mappingType='discrete', mappingColumn='border_width', discrete=[str_width, str_width])
        mappings.append(border_width_map)

    if 'border_color' in statistics:
        border_color = statistics['border_color']
        str_color = [str(x) for x in np.unique(border_color)]
        border_color_map = cytoscape.vizmap.mapVisualProperty(
            visualProperty='NODE_BORDER_PAINT', mappingType='discrete', mappingColumn='border_color', discrete=[str_color, str_color])
        mappings.append(border_color_map)

    if 'impact_factor' in statistics:
        impact_factor = list(statistics['impact_factor'].values())
        # node_transparency = cytoscape.vizmap.mapVisualProperty(
        #    visualProperty='NODE_TRANSPARENCY', table='node', mappingType='continuous', mappingColumn='impact_factor', lower=[min(impact_factor), 10], center=[np.mean(impact_factor), 135], upper=[max(impact_factor), 255])
        # mappings.append(node_transparency)

    node_size = cytoscape.vizmap.mapVisualProperty(
        visualProperty='NODE_SIZE', table='node', mappingType='continuous', mappingColumn='weight', lower=[1, 17], center=[5, 30], upper=[100, 80])
    mappings.append(node_size)

    update_cytoscape_style(style_name, mappings)
    time.sleep(10)


def convert_genes(genes, conv_from, conv_to, species):
    mg = get_client('gene')
    out = {}
    i = 0
    tmp_fname = 'data/enst-tmp.pickle'
    if os.path.isfile(tmp_fname):
        out = pickle.load(open(tmp_fname, 'rb'))
    for gene in genes:
        i += 1
        if gene in out:
            continue
        fields = ','.join([conv_to, 'refseq'])
        query = mg.query(gene, scopes=conv_from,
                         fields=fields, species=species)
        out[gene] = None
        if 'hits' in query:
            if len(query['hits']) > 0:
                if 'symbol' in query['hits'][0]:
                    out[gene] = query['hits'][0]['symbol']
        if i % 1000 == 0:
            pickle.dump(out, open(tmp_fname, 'wb'))
            print(i, gene, out[gene])
    pickle.dump(out, open(tmp_fname, 'wb'))
    mapping = out
    no_results = []
    new_ids = []
    for x in genes:
        new_ids.append(mapping[x])
        if mapping[x] is None:
            no_results.append(x)
    # for x in out:
    #    if conv_to in x.keys():
    #        mapping[x['query']] = x[conv_to]
    #    else:
    #        no_results.append(x)
    # for x in genes:
    #    if x not in mapping.keys():
    #        no_results.append(x)
    # new_ids = [mapping[x] if x not in no_results else '-' for x in genes]
    # if '-' in new_ids:
    #    new_ids.remove('-')
    return new_ids, mapping, no_results


def run_gene_scf(genes):
    geneSCF_folder = 'geneSCF'
    os.chdir(geneSCF_folder)
    data_folder = 'data'
    gene_list_file = 'gene_list.csv'
    genes = pd.DataFrame(genes)
    # genes = pd.DataFramed([x.upper() for x in genes])
    genes.to_csv('{}/{}'.format(data_folder, gene_list_file),
                 header=None, index=None)
    os.system(
        './geneSCF -m=normal -i={}/{} -t=sym -o={}/out -db=KEGG -bg=20000 -p=no -org=mmu'.format(data_folder, gene_list_file, data_folder))


def weigh_edges(G):
    weight_map = {-1: 0.05, 0: -0.1, 1: 0.4}
    weight_map = {-1: 0.3, 0: 0.3, 1: 0.3}
    for node1, node2, idx in G.edges:
        if 'predictions' not in G[node1][node2][idx]:
            break
        G[node1][node2][idx]['weight'] = weight_map[G[node1]
                                                    [node2][idx]['predictions']]
    return G


def edge_weight(G, node1, node2, base=2, max_val=100):
    # if G.node[node1]['weight'] < 2 or G.node[node2]['weight'] < 2:
    degree1 = len(G.node[node1])
    degree2 = len(G.node[node2])
    return float(degree1*degree2/max([np.power(degree1-degree2, 2), 1]))/max_val  # float(min(10./np.power(base, max([G.node[node1]['weight'], G.node[node2]['weight']])), max_val))/max_val
    # return float(min(10./np.power(base, max([G.node[node1]['weight'], G.node[node2]['weight']])), max_val))/max_val


def get_communities_grivan_newman(G):
    communities_generator = community.girvan_newman(G)
    next_level_communities = None
    for comm in communities_generator:
        max_clique = max([len(x) for x in comm])
        print(max_clique)
        if max_clique < 0.1 * len(G.node) or max_clique < 15:
            next_level_communities = list(comm)
            break
    return next_level_communities


def get_communities_greedy(G):
    communities_generator = community.greedy_modularity_communities(G)
    next_level_communities = list(communities_generator)
    return next_level_communities


def get_communities_k_clique(G, k):
    return list(community.k_clique_communities(G, k))


def get_communities_async_label_propagation(G):
    return list(community.asyn_lpa_communities(G))


def get_communities_label_propagation(G):
    return list(community.label_propagation_communities(G))


def get_communities_fluid(G):
    connected_components = components.connected_components(G)
    modules = []
    min_size = 50
    coef = 1./min_size
    for component in connected_components:
        if len(component) < min_size:
            modules = modules + [component]
            continue
        k = int(np.ceil(coef * len(G.node)))
        modules = modules + list(community.asyn_fluidc(G.subgraph(component), k, seed=123))
    return modules


def weight_by_community(G):
    def is_big(G, node):
        return len(G.node[node]) > 10
    # next_level_communities = get_communities_grivan_newman(G)
    # next_level_communities = get_communities_greedy(G)
    # next_level_communities = get_communities_k_clique(G, 4)
    # next_level_communities = get_communities_label_propagation(G)
    next_level_communities = get_communities_fluid(G)
    # print(next_level_communities)
    for i in range(len(next_level_communities)-1):
        community1 = next_level_communities[i]
        for j in np.arange(i+1, len(next_level_communities)):
            community2 = next_level_communities[j]

            for node1 in community1:
                if not is_big(G, node1):
                    continue
                for node2 in community2:
                    if node2 not in list(G[node1].keys()):
                        if not is_big(G, node2):
                            continue
                        G.add_edge(node1, node2)
                        G[node1][node2][0]['betweenness'] = 0
                        G[node1][node2][0]['transparency'] = 0.
                        G[node1][node2][0]['predictions'] = -100
                        G[node1][node2][0]['weight'] = 0.1  # edge_weight(G, node1, node2, 2, 100)
    return G


def weight_by_degree(G):
    def is_big(G, node1):
        return len(G[node1]) > 3

    for i in range(len(nodes)-1):
        node1 = nodes[i]
        if not is_big(G, node1):
            continue
        for j in np.arange(i+1, len(G.nodes)):
            node2 = nodes[j]
            if not is_big(G, node2):
                continue
            if node2 not in list(G[node1].keys()):
                G.add_edge(node1, node2)
                G[node1][node2][0]['betweenness'] = 0
                G[node1][node2][0]['transparency'] = 0.
                G[node1][node2][0]['predictions'] = -100
                G[node1][node2][0]['weight'] = edge_weight(G, node1, node2, 2, 100)
    return G


def weigh_module(G):
    for node in G.nodes:
        G.node[node]['weight'] = len(G[node])
    for node1, node2, idx in G.edges:
        G[node1][node2][idx]['weight'] = edge_weight(G, node1, node2, 2, 100)
    nodes = list(G.nodes.keys())

    # G = weight_by_degree(G)
    G = weight_by_community(G)

    return G


def pathway_map(G, reverse_mapping):
    pathways = pickle.load(open('pathways.pickle', 'rb'))
    for pathway_name, pathway_genes in pathways.items():
        for gene in pathway_genes:
            if gene is not None and gene in G.nodes:
                if 'pathway' in G.node[gene]:
                    G.node[gene]['pathway'] = '{} ; {}'.format(
                        G.node[gene]['pathway'], pathway_name)
                else:
                    G.node[gene]['pathway'] = pathway_name
    return G


def calculate_edge_betweenness_centrality(G):
    filename = 'data/edge_betweenness.pickle'
    betweenness = None
    if not os.path.isfile(filename):
        betweenness = nx.edge_betweenness_centrality(G, k=len(G.nodes))
        pickle.dump(betweenness, open(filename, 'wb'))
    else:
        betweenness = pickle.load(open(filename, 'rb'))
    return betweenness


def map_betweenness_centrality_to_subnetwork(G_from, G_to):
    betweenness = calculate_edge_betweenness_centrality(G_from)
    for node1, node2, idx in G_to.edges:
        key = None
        if (node1, node2) in betweenness:
            key = (node1, node2)
        elif (node2, node1) in betweenness:
            key = (node2, node1)
        if key is None:
            continue
        G_to[node1][node2][idx]['betweenness'] = betweenness[key]
        G_to[node2][node1][idx]['betweenness'] = betweenness[key]

    return G_to, betweenness


def find_hubs(G, degree_cutoff=70):
    hubs = []
    degrees = {node: len(G[node]) for node in G.nodes}
    for node, degree in degrees.items():
        if degree > degree_cutoff:
            hubs.append(node)
    # sorted_degrees = sorted(degrees.values(), reverse=True)
    # plt.plot(sorted_degrees)
    # plt.show()
    return set(hubs)


def map_hubs_to_subnetwork(G_from, G_to):
    hubs = find_hubs(G_from)
    for node in G_to.nodes:
        if node in hubs:
            G_to.node[node]['hub'] = True
        else:
            G_to.node[node]['hub'] = False
    return G_to


def add_node_shape(G):
    for node in G.nodes:
        G.node[node]['shape'] = 'circle'
        if 'color' not in G.node[node]:
            continue
        if G.node[node]['color'] == 'red' or G.node[node]['color'] == 'yellow':
            G.node[node]['shape'] = 'triangle'
        if 'hub' in G.node[node] and (G.node[node]['color'] == 'red' or G.node[node]['color'] == 'yellow'):
            if G.node[node]['hub']:
                G.node[node]['shape'] = 'hexagon'
    return G


def add_transparency(G):
    low = 50
    high = 255
    addition = 0
    betweenness = []
    for node1, node2, idx in G.edges:
        betweenness.append(G[node1][node2][idx]['betweenness'])

    bet_min = min(betweenness)
    bet_max = max(betweenness)

    transparency = []
    for node1, node2, idx in G.edges:
        G[node1][node2][idx]['transparency'] = int(
            (G[node1][node2][idx]['betweenness']-bet_min)/(bet_max - bet_min)*(high-low)+low)
        if 'predictions' in G[node1][node2][idx]:
            if G[node1][node2][idx]['predictions'] in [-1, 1]:
                G[node1][node2][idx]['transparency'] = min(
                    [high, G[node1][node2][idx]['transparency']+addition])
        transparency.append(G[node1][node2][idx]['transparency'])
    return G, transparency


def add_node_boundary(G, pathway_border_colors):
    pathway_border_width = 2
    hub_border_width = 5
    hub_border_color = '#3d2aeb'

    border_widths = [pathway_border_width, hub_border_width]
    border_colors = [hub_border_color]
    for node in G.nodes:
        if 'pathway' in G.node[node]:
            pass
            # G.node[node]['border_width'] = pathway_border_width
            # pathway_color = pathway_border_colors[G.node[node]['pathway']]
            # G.node[node]['border_color'] = pathway_color
            # border_colors.append(pathway_color)
        if G.node[node]['hub']:
            G.node[node]['border_width'] = hub_border_width
            G.node[node]['border_color'] = hub_border_color

    return G, border_widths, border_colors


def get_pathway_border_colors():
    colors = {}
    colors['t2d'] = '#aa48db'
    return colors


def add_node_names(G, names=None):
    for node in G.nodes:
        if names is None:
            G.node[node]['name'] = node
        else:
            G.node[node]['name'] = names[node]
    return G


def add_impact_factor(G, impact_factor):
    for node in G.nodes:
        if node not in impact_factor:
            G.node[node]['impact_factor'] = 0.0
            # print(node)
        else:
            G.node[node]['impact_factor'] = impact_factor[node]
    return G


def preprocess_network(G, G_mapping_from=None):
    # new_ids, mapping, no_results = convert_genes(
    #    list(G_diabetes_full.nodes), 'symbol', 'entrezgene', 'mouse')
    # reverse_mapping = {v: k for k, v in mapping.items()}

    statistics = {}
    if G_mapping_from is not None:
        print(G_full is None)
        G = map_hubs_to_subnetwork(G_mapping_from, G)
        G, betweenness = map_betweenness_centrality_to_subnetwork(
            G_mapping_from, G)
        G, transparency = add_transparency(G)
        statistics['betweenness'] = betweenness
        statistics['transparency'] = transparency

    pathway_border_colors = get_pathway_border_colors()
    G = add_node_shape(G)
    # G = prune_network(G, [], [], [
    #    'predictions'], [[-1, 1]], hub_threshold=15, affect_neighbors=False)
    # G = pathway_map(G, reverse_mapping)
    G, border_width, border_colors = add_node_boundary(
        G, pathway_border_colors)
    statistics['border_width'] = border_width
    statistics['border_color'] = border_colors
    G = weigh_edges(G)
    G = add_node_names(G)
    return G, statistics


def load_impact_factor_gene_iso(bs, imp):
    isoforms = bs['isoform']
    groups = set([x.split('-')[0] for x in isoforms])
    impact_factor_file = "data/mock_IF.pickle"

    impact_factor = None
    if not os.path.isfile(impact_factor_file):
        max_size = -1
        for group in groups:
            size = len([i for i, s in enumerate(isoforms) if group == s.split('-')[0]])
            if size > max_size:
                max_size = size

        print("Max size: {}".format(max_size))
        impact_factor = {}
        for group in groups:
            idx = [i for i, s in enumerate(isoforms) if group == s.split('-')[0]]
            to_pad = max_size - len(idx)
            # bs_change = np.concatenate(bs['bs'][idx].values)
            # print(len(bs['bs'][idx].values))
            if len(bs['bs'][idx].values) > 100:
                print(isoforms[idx], len(bs['bs'][idx].values), group)
            impact_factor[group] = [imp.gene_if(1, -1, -1, x) for x in bs['bs'][idx].values] + list(-1*np.ones(to_pad))
            # print(impact_factor[group])
        pickle.dump(impact_factor, open(impact_factor_file, 'wb'))
    else:
        impact_factor = pickle.load(open(impact_factor_file, 'rb'))
    return impact_factor


def load_impact_factor_gene_tissue(bs, imp):
    isoforms = bs['isoform']
    print(isoforms)
    impact_factor_file = "data/mock_IF_tissue.pickle"
    expression = load_tissue_level()
    expression = expression.groupby(['symbol']).sum()
    expression = expression.drop([0])
    tissues = expression.columns
    genes = set([x.split('-')[0] for x in isoforms])

    impact_factor = None
    if not os.path.isfile(impact_factor_file):
        i = 0
        impact_factor = {}
        for gene in genes:
            if gene not in expression.index:
                print(gene)
                i += 1
                continue
            idx = [i for i, s in enumerate(isoforms) if gene == s.split('-')[0]]

            bs_change = np.concatenate(bs['bs'][idx].values)
            max_expr = max(expression.loc[gene])
            impact_factor[gene] = [imp.gene_if(expression.loc[gene, tissue]/max_expr, -1, -1, bs_change) for tissue in tissues]
            # print(impact_factor[group])
        pickle.dump(impact_factor, open(impact_factor_file, 'wb'))
        print('Not found genes: {}'.format(i))
    else:
        impact_factor = pickle.load(open(impact_factor_file, 'rb'))
    return impact_factor


def map_compas_ids():
    # from Bio import SeqIO
    mapping_file = "data/enst-compas-mapping.pickle"

    compas_file = "data/AS_DB_Final.fa"
    biomart_file = "data/enst-biomart.txt"

    compas_enst_mapping = {}
    if os.path.isfile(mapping_file):
        compas_enst_mapping = pickle.load(open(mapping_file, 'rb'))
        return compas_enst_mapping

    compas_ids = pd.read_csv("data/isoform_seq_uniq_id.csv", header=None, sep=',', index_col=None)
    compas_ids_dict_split = compas_ids.to_dict('records')
    compas_ids_dict = {}
    count = 0
    for record in compas_ids_dict_split:
        if 'ENST' in record[2]:
            enst = [x for x in record[2].replace('|', ',').split(',') if 'ENST' in x][0].split('.')[0].split(' ')[0].split('\t')[0].rstrip()
            compas_enst_mapping[record[0]] = enst
            count += 1
        compas_ids_dict[record[1]] = record[0]
    print('ENST in COMP-AS: {}'.format(count))

    biomart_dict = {}
    for record in SeqIO.parse(biomart_file, "fasta"):
        enst = record.id.split('|')[2]
        biomart_dict[str(record.seq[:-1])] = enst

    i = 0
    j = 0
    k = 0
    compas_seq_dict = {}
    for record in SeqIO.parse(compas_file, "fasta"):
        seq = str(record.seq)
        if 'ENST' not in record.id:
            if seq in biomart_dict:
                compas_seq_dict[seq] = biomart_dict[seq]
                j += 1
            else:
                k += 1
        else:
            enst = [x for x in record.id.replace('|', ',').split(',') if 'ENST' in x][0].split('.')[0]
            i += 1
            compas_seq_dict[seq] = enst

    print('Compas records with ENST: {}'.format(i))
    print('Biomart records: {}'.format(j))
    print('Unknown records: {}'.format(k))

    from Bio import pairwise2

    identity_threshold = 0.9

    filtered_file = 'data/filtered_seq.pickle'
    filtered_out = {}
    if os.path.isfile(filtered_file):
        filtered_out = pickle.load(open(filtered_file, 'rb'))
    counter = 10

    for seq in compas_ids_dict.keys():

        continue

        if compas_ids_dict[seq] in compas_enst_mapping:
            continue

        if seq in filtered_out:
            continue

        if seq in compas_seq_dict:
            compas_enst_mapping[compas_ids_dict[seq]] = compas_seq_dict[seq]
            counter += 1
            continue

        from Bio.Blast.Applications import NcbiblastpCommandline
        from io import StringIO
        from Bio.Blast import NCBIXML
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        best_seq = None
        best_identity = -1

        seq1 = SeqRecord(Seq(seq), id="seq1")
        SeqIO.write(seq1, "data/seq1.fasta", "fasta")
        output = NcbiblastpCommandline(query="data/seq1.fasta", subject=biomart_file, outfmt=5)()[0]
        blast_result_record = NCBIXML.read(StringIO(output))

        print(blast_result_record)
        print(blast_result_record.alignments[0].title)

        # for alignment in blast_result_record.alignments:
        #    for hsp in alignment.hsps:
        #        print(hsp.expect)
        #        break
        # print(hsp.identities/len(hsp.sbjct.strip('-')))
        # print(hsp.sbjct)

        def check_template():
            for template_seq in biomart_dict.keys():
                n = len(template_seq)
                m = len(seq)
                if np.min((n, m))/max((n, m)) < identity_threshold:
                    continue

                for alignment in pairwise2.align.globalxs(template_seq, seq, -10, -0.5):
                    matches = sum(aa1 == aa2 for aa1, aa2 in zip(alignment[0], alignment[1]))
                    identity = float(matches) / max((n, m))
                    if identity > best_identity:
                        best_seq = template_seq
                        best_identity = identity

            counter += 1
            print(counter, best_identity)
            if best_identity >= identity_threshold:
                print(counter, best_identity)
                compas_enst_mapping[compas_ids_dict[seq]] = biomart_dict[best_seq]
            else:
                filtered_out[seq] = '{}-{}'.format(best_identity, biomart_dict[best_seq])

        if counter % 10 == 0:
            print(counter)
            pickle.dump(compas_enst_mapping, open(mapping_file, 'wb'))
            pickle.dump(filtered_out, open(filtered_file, 'wb'))
    pickle.dump(compas_enst_mapping, open(mapping_file, 'wb'))

    print(len(compas_enst_mapping))
    return compas_enst_mapping


def load_if_color_scale(impact_factor):
    orig_impact_factor = list(impact_factor.values())
    impact_factor = sorted(orig_impact_factor)
    n = len(orig_impact_factor)
    top_5_val = impact_factor[-int(0.05*n)]
    top_15_val = impact_factor[-int(0.15*n)]
    top_25_val = impact_factor[-int(0.25*n)]
    top_50_val = impact_factor[-int(0.5*n)]
    top_75_val = impact_factor[-int(0.75*n)]
    top_100_val = impact_factor[0]
    vals = [top_5_val, top_15_val, top_25_val, top_50_val, top_75_val, top_100_val]
    colors = ['#420D09', '#800000', '#FA8072', '#1890F0', '#A8D8F0', '#ebf4f5']
    return vals, colors


def plot_clustermap(df_matrix, fname, colors):
    from bokehheat import heat, jheat
    mx = df_matrix.max().max()
    o_clustermap, ls_xaxis, ls_yaxis = heat.clustermap(
        df_matrix=df_matrix,
        ls_color_palette=colors,
        r_low=0,
        r_high=mx,
        s_z="log2",
        # tt_axis_annot=tt_boolecatquant,
        b_ydendo=True,
        b_xdendo=True,
        # s_method='average',
        # s_metric='euclidean',
        # b_optimal_ordering=True,
        #i_px = 64,
        #i_height = 12,
        #i_width = 12,
        #i_min_border_px = 128,
        s_filename=fname,
        s_filetitel="the Clustermap",
    )
    save(o_clustermap)


def impact_factor_profile(G, out_file):
    factors = G.nodes

    impact_factor = load_impact_factor()  # load_impact_factor_gene_tissue(bs, imp)
    genes_if = impact_factor['genes_max']
    impact_factor = impact_factor['genes_tissue_max']
    impact_factor = impact_factor.drop([x for x in impact_factor.index if x not in factors])
    if impact_factor.empty:
        return
    if_all = impact_factor.values.flatten()
    print(impact_factor)
    # impact_factor = pd.DataFrame.from_dict(impact_factor).transpose()
    # print(impact_factor)

    tissues = list(impact_factor.columns)
    genes = list(impact_factor.index)
    # gtex = pd.read_csv('data/tissues_SMTS_tpm_mean_curated.tsv', sep=',')
    # tissues = list(gtex.columns[1:].values)
    # impact_factor.columns = tissues
    # print(tissues)
    # tissues = [list(x.keys())[0] for x in impact_factor.loc[genes[0]]]
    # impact_factor.columns = tissues
    # tissues = list(np.array(impact_factor.columns, dtype=str))

    df = pd.DataFrame(impact_factor.stack()).reset_index()
    print(df)
    df.columns = ['Gene', 'Tissue', 'IF']
    # colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce", "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]
    # mapper = LinearColorMapper(palette=colors, low=np.min(if_all), high=np.max(if_all))  # df.IF.min()

    TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"

    p = figure(title="Impact Factor Profile",
               x_range=tissues, y_range=genes,
               x_axis_location="above", plot_width=400, plot_height=900,
               tools=TOOLS, toolbar_location='below',
               tooltips=[('Impact factor', '@Tissue @Gene'), ('IF', '@IF'), ('Percentile', '@Percentile')])
    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "7px"
    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = np.pi / 3

    fill_colors = []
    percentiles = []
    thresholds, colors = load_if_color_scale(genes_if)
    threshold_percentile = ['5%', '15%', '25%', '50%', '75%', '100%']
    for impact_factor_value in df.IF:
        for threshold, col, perc in zip(thresholds, colors, threshold_percentile):
            if impact_factor_value >= threshold:
                fill_colors.append(col)
                percentiles.append(perc)
                break
    df['Color'] = fill_colors
    df['Percentile'] = percentiles
    p.rect(x="Tissue", y="Gene", width=1, height=1,
           source=df,
           fill_color="Color",  # {'field': 'IF', 'transform': } fill_colors,
           line_color=None)
    # {'field': 'IF', 'transform': mapper},
    # color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="7px",
    #                     ticker=BasicTicker(desired_num_ticks=len(colors)),
    #                     formatter=PrintfTickFormatter(format="%d%%"),
    #                     label_standoff=6, border_line_color=None, location=(0, 0))
    #p.add_layout(color_bar, 'right')
    output_file(out_file)
    save(p)
    plot_clustermap(impact_factor, '{}-cluster.html'.format(out_file), colors[::-1])


def color_by_if(G, impact_factor):
    vals, colors = load_if_color_scale(impact_factor)
    for node in G.node:
        G.node[node]['color'] = colors[-1]
        for val, col in zip(vals, colors):
            if G.node[node]['impact_factor'] >= val:
                G.node[node]['color'] = col
                break
    return G


if __name__ == "__main__":
    graph_file = 'data/graph.pickle'
    statistics_file = 'data/statistics.pickle'
    modules = load_disease_modules()
    impact_factor = load_impact_factor()

    G_full = None
    statistics = None
    if not os.path.isfile(graph_file):
        interactome = load_interactome()

        cols = interactome.columns
        G_full = nx.from_pandas_edgelist(
            interactome[cols], source=cols[0], target=cols[1], create_using=nx.MultiGraph)

        loops = list(G_full.selfloop_edges())
        G_full.remove_edges_from(loops)
        for node1, node2, idx in G_full.edges:
            G_full[node1][node2][idx]['predictions'] = 1
        # print(G_full is None)
        G_full, statistics = preprocess_network(G_full, G_full)
        pickle.dump(G_full, open(graph_file, 'wb'))
        pickle.dump(statistics, open(statistics_file, 'wb'))
    else:
        G_full = pickle.load(open(graph_file, 'rb'))
        statistics = pickle.load(open(statistics_file, 'rb'))

    # print(statistics)
    statistics['impact_factor'] = impact_factor['genes_max']
    # statistics['impact_factor_color'] =
    G_full = add_impact_factor(G_full, impact_factor['genes_max'])
    cytoscape = cyrest.cyclient()

    for key in modules.keys():
        if 'Dermatitis, Occupational' not in key:
            continue
        out_fname = 'disease-figures/{}.pdf'.format(key)
        print(out_fname)
        # if os.path.isfile(out_fname):
        #    continue
        print(modules[key])
        print(G_full.nodes)
        G_module = prune_network(deepcopy(G_full), ['name'], [modules[key]], [], [], hub_threshold=10000, affect_neighbors=True)
        G_module = color_nodes(G_module, modules[key], 'red')
        G_module = add_node_shape(G_module)
        G_module = color_by_if(G_module, impact_factor['genes_max'])
        G_module = weigh_module(G_module)

        break
        continue
        impact_factor_profile(G_module, "disease-profiles/{}.html".format(key))

        graph = nx.readwrite.json_graph.cytoscape_data(G_module)

        cytoscape.session.new()
        network_id, res = create_new_network(graph)
        apply_diabetes_style(cytoscape, statistics)
        time.sleep(2)
        # cytoscape.layout.force_directed(EdgeAttribute='weight', maxWeightCutoff=100.0, defaultNodeMass=2, defaultSpringCoefficient=5e-5,
        #                                Type='normalized value', numIterations='1000')
        cytoscape.layout.apply_preferred()
        res = cytoscape.result(filetype="PDF", saveas=out_fname)
        time.sleep(2)
        export_fname = os.path.join(os.path.abspath('.'), 'disease-networks/{}.{}')
        # export_fname = 'disease-networks/{}.{}'
        # cytoscape.network.export(options='CYJS', OutputFile=export_fname.format(key, 'js'))
        # cytoscape.vizmap.export(options='json', OutputFile=export_fname.format('cytoscape-style', 'json'))

        # cytoscape.network.deselect(nodeList='all', edgeList='all')
        # cytoscape.view.fit_content()
        # cytoscape.view.export(options="PDF", outputFile=os.path.join(os.path.abspath('.'), 'disease-figures/{}.pdf'.format(key)))
        # cytoscape.result(filetype="CYS", saveas=export_fname)
        # cytoscape.network.export(OutputFile=export_fname, verbose=True)
        # cytoscape.network.export(OutputFile=export_fname, options='CX JSON (.cx)', verbose=True)
        print(res)
        time.sleep(10)
        res = cytoscape.session.runGarbageCollection()
        print(res)
        time.sleep(10)
