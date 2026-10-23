# -*- coding: utf-8 -*-
"""抓取正式基因清单：UniProt分泌蛋白 + KEGG免疫通路 + KEGG PANoptosis"""
import urllib.request, json, time, re, os

DATA = r'C:\Users\98757\Desktop\生信分析\data'

def get(url, binary=False):
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    d = urllib.request.urlopen(req, timeout=60).read()
    return d if binary else d.decode('utf-8', 'ignore')

def kegg_genes(pid):
    txt = get(f'https://rest.kegg.jp/get/{pid}')
    genes = []
    in_gene = False
    for line in txt.splitlines():
        if line.startswith('GENE'):
            in_gene = True
            line = line[12:]
        elif in_gene and line.startswith('            '):
            pass
        else:
            in_gene = False
            continue
        m = re.match(r'\s*\d+\s+([A-Z0-9]+)[;(\s]', line)
        if m:
            genes.append(m.group(1))
    return genes

# 1. KEGG 免疫通路合集（选题B）
IMMUNE_PATHWAYS = ['hsa04060','hsa04061','hsa04062','hsa04610','hsa04620','hsa04621','hsa04622',
                   'hsa04623','hsa04625','hsa04612','hsa04640','hsa04650','hsa04657','hsa04658',
                   'hsa04659','hsa04660','hsa04662','hsa04664','hsa04666','hsa04670','hsa04672']
immune = set()
for p in IMMUNE_PATHWAYS:
    g = kegg_genes(p)
    print(p, len(g))
    immune.update(g)
    time.sleep(0.3)
print('免疫基因合计:', len(immune))

# 2. KEGG PANoptosis 合集（选题C）：坏死性凋亡+凋亡+NOD样受体（焦亡机器）
PANOPT_PATHWAYS = ['hsa04217', 'hsa04210', 'hsa04621']
panopt = set()
for p in PANOPT_PATHWAYS:
    g = kegg_genes(p)
    print(p, len(g))
    panopt.update(g)
    time.sleep(0.3)
# 补充 PANoptosome 关键分子（KEGG 通路未完整收录）
panopt.update(['ZBP1','AIM2','NLRP3','NLRC4','PYCARD','GSDMD','GSDME','CASP1','CASP6','MLKL',
               'RIPK1','RIPK3','TNF','FAS','FASLG','CYLD','TNFRSF1A','TRADD','FADD'])
print('PANoptosis 基因合计:', len(panopt))

# 先保存已成功的清单
json.dump({'immune_kegg': sorted(immune), 'panoptosis_kegg': sorted(panopt)},
          open(os.path.join(DATA, 'gene_lists.json'), 'w', encoding='utf-8'), indent=1)
print('KEGG 清单已保存')

# 3. UniProt 分泌蛋白（KW-0964, human, reviewed），断流自动重试
import http.client
url = ('https://rest.uniprot.org/uniprotkb/stream?query=keyword:KW-0964%20AND%20organism_id:9606'
       '%20AND%20reviewed:true&format=tsv&fields=accession,gene_names')
tsv = None
for attempt in range(5):
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        tsv = urllib.request.urlopen(req, timeout=90).read().decode('utf-8', 'ignore')
        break
    except http.client.IncompleteRead as e:
        tsv = e.partial.decode('utf-8', 'ignore')
        if tsv.count('\n') > 3000:
            break
        time.sleep(3)
    except Exception as e:
        print('uniprot err', e); time.sleep(3)
secreted = set()
for i, line in enumerate(tsv.splitlines()):
    if i == 0 or not line.strip():
        continue
    parts = line.split('\t')
    if len(parts) >= 2 and parts[1].strip():
        secreted.add(parts[1].split()[0])
print('UniProt 分泌蛋白合计:', len(secreted))

json.dump({'immune_kegg': sorted(immune), 'panoptosis_kegg': sorted(panopt), 'secreted_uniprot': sorted(secreted)},
          open(os.path.join(DATA, 'gene_lists.json'), 'w', encoding='utf-8'), indent=1)
print('已保存 gene_lists.json')
