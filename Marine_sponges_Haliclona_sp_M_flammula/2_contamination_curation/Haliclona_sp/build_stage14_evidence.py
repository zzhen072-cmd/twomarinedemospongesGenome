#!/usr/bin/env python3import csv, json, math, statisticsfrom pathlib import PathROOT = Path('/work1/sponge/lht_reassembly_compare')ADAPTOR = ROOT/'10_fcs_gx/lht_fcs_adaptor.fasta'GX = ROOT/'10_fcs_gx/lht_fcs_gx.fasta'TIARA = ROOT/'13_tiara/lht.tiara.tsv'KRAKEN = ROOT/'11_kraken2/lht.kraken2_c0.out'BLOB = ROOT/'12_blobtoolkit/lht_fcs_adaptor.blobdir'DB = Path('/public/home/sponge/sci_data/S07/software/standard16')OUT = ROOT/'14_manual_decontam'OUT.mkdir(exist_ok=True)def norm(s): return s.strip().split()[0].lstrip('>')def fasta_records(path):    name, header, seq = None, None, []    with open(path) as fh:        for line in fh:            if line.startswith('>'):                if name is not None: yield name, header, ''.join(seq)                header = line[1:].rstrip(); name = norm(header); seq = []            else: seq.append(line.strip())    if name is not None: yield name, header, ''.join(seq)records = {n:(h,s) for n,h,s in fasta_records(ADAPTOR)}gx_ids = {n for n,_,_ in fasta_records(GX)}

tiara = {}
with open(TIARA) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        tiara[norm(row['sequence_id'])] = row

parent, rank, taxname = {}, {}, {}
taxdir = DB
with open(taxdir/'nodes.dmp') as fh:
    for line in fh:
        p=line.split('|'); parent[int(p[0])] = int(p[1]); rank[int(p[0])] = p[2].strip()
with open(taxdir/'names.dmp') as fh:
    for line in fh:
        p=line.split('|')
        if p[3].strip() == 'scientific name': taxname[int(p[0])] = p[1].strip()

def domain(t):
    seen=set()
    while t and t not in seen:
        seen.add(t)
        if t == 2: return 'bacteria'
        if t == 2157: return 'archaea'
        if t == 2759: return 'eukarya'
        if t == 10239: return 'viruses'
        q=parent.get(t,t)
        if q == t: break
        t=q
    return 'other'

kraken = {}
with open(KRAKEN) as fh:
    for line in fh:
        p=line.rstrip('\n').split('\t')
        if len(p) >= 4:
            tid = int(p[2]) if p[2].isdigit() else 0
            kraken[norm(p[1])] = (p[0], tid, taxname.get(tid,''), domain(tid))

def vector(name):
    x=json.load(open(BLOB/name))
    if isinstance(x,dict):
        for key in ('values','data'):
            if key in x: return x[key]
    if isinstance(x,list): return x
    raise ValueError(f'Unsupported BlobToolKit JSON: {name}: {list(x)}')

ids=[norm(str(x)) for x in vector('identifiers.json')]
length=vector('length.json'); gc=vector('gc.json'); cov=vector('lht.ont_to_fcs_adaptor_cov.json')
blob={ids[i]:(int(length[i]),float(gc[i]),float(cov[i])) for i in range(len(ids))}

host_cov=[v[2] for k,v in blob.items() if tiara.get(k,{}).get('class_fst_stage')=='eukarya' and k in gx_ids and v[2]>0]
host_gc=[v[1] for k,v in blob.items() if tiara.get(k,{}).get('class_fst_stage')=='eukarya' and k in gx_ids]
med_cov=statistics.median(host_cov); med_gc=statistics.median(host_gc)
mad_cov=statistics.median(abs(math.log2(x/med_cov)) for x in host_cov)
mad_gc=statistics.median(abs(x-med_gc) for x in host_gc)

rows=[]; delete=[]; review=[]
for name,(header,seq) in records.items():
    L,G,C=blob.get(name,(len(seq),100*(seq.upper().count('G')+seq.upper().count('C'))/max(1,len(seq)),0))
    tr=tiara.get(name,{})
    t1=tr.get('class_fst_stage','missing'); t2=tr.get('class_snd_stage','n/a')
    ks,kt,kn,kd=kraken.get(name,('U',0,'','unclassified'))
    fcs=name not in gx_ids
    pro_t=t1 in {'bacteria','archaea','prokarya'}
    pro_k=kd in {'bacteria','archaea','viruses'}
    cov_out=C>0 and abs(math.log2(C/med_cov)) > max(2.5, 6*mad_cov)
    gc_out=abs(G-med_gc) > max(0.12 if med_gc <= 1 else 12, 6*mad_gc)
    support=sum((fcs,pro_t,pro_k,cov_out,gc_out))
    conflict=(t1=='eukarya' and kd=='eukarya') or (t1=='eukarya' and not pro_k)
    if (fcs and (pro_t or pro_k or (cov_out and gc_out))) or (pro_t and pro_k) or (pro_t and cov_out and gc_out) or (pro_t and gc_out and kd != 'eukarya') or (t1 == 'eukarya' and pro_k and gc_out):
        decision='DELETE' if not (fcs and conflict and support<3) else 'REVIEW'
    elif fcs or pro_t or pro_k or cov_out or gc_out:
        decision='REVIEW'
    else: decision='KEEP'
    reason=';'.join(k for k,v in [('fcs_gx_removed',fcs),('tiara_prok',pro_t),('kraken_prok_or_virus',pro_k),('coverage_outlier',cov_out),('gc_outlier',gc_out)] if v) or 'no_contamination_signal'
    if decision=='DELETE': delete.append(name)
    elif decision=='REVIEW': review.append(name)
    rows.append([name,L,G,C,fcs,t1,t2,ks,kt,kn,kd,cov_out,gc_out,support,decision,reason])

head=['contig','length','gc','ont_cov','fcs_gx_removed','tiara_stage1','tiara_stage2','kraken_status','kraken_taxid','kraken_name','kraken_domain','coverage_outlier','gc_outlier','support_count','decision','reason']
with open(OUT/'lht_decontam_evidence.tsv','w',newline='') as fh:
    w=csv.writer(fh,delimiter='\t'); w.writerow(head); w.writerows(rows)
for fn,items in [('delete_candidates.txt',delete),('review_candidates.txt',review)]:
    (OUT/fn).write_text(''.join(x+'\n' for x in items))
with open(OUT/'audit_summary.txt','w') as fh:
    fh.write(f'input_contigs\t{len(records)}\ninput_bases\t{sum(len(s) for _,s in records.values())}\n')
    fh.write(f'fcs_gx_removed\t{len(set(records)-gx_ids)}\ndelete_candidates\t{len(delete)}\nreview_candidates\t{len(review)}\n')
    fh.write(f'delete_bases\t{sum(len(records[x][1]) for x in delete)}\nreview_bases\t{sum(len(records[x][1]) for x in review)}\n')
    fh.write(f'host_median_cov\t{med_cov}\nhost_median_gc\t{med_gc}\n')
print((OUT/'audit_summary.txt').read_text(), end='')

delete_set=set(delete)
with open(OUT/'lht_manual_decontaminated.fasta','w') as keep, open(OUT/'lht_removed_contaminants.fasta','w') as removed:
    for name,(header,seq) in records.items():
        fh = removed if name in delete_set else keep
        fh.write('>'+header+'\n')
        for i in range(0,len(seq),80): fh.write(seq[i:i+80]+'\n')
