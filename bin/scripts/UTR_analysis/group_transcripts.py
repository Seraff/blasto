#!/usr/bin/env python3

from fileinput import filename
from os.path import exists
import json

from pprint import pprint
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import intervals as I

INPUT_FOLDER = 'input_06_07'
TRANSCRIPT_FILENAME = f'{INPUT_FOLDER}/transcript.fasta'
BLAST_FILENAME = f'{INPUT_FOLDER}/blast_to_genome.tsv'

GENOME_FILENAME = 'genome_good.fasta'
ANNOTATION_FILENAME = 'annotation.gff'

# Output
OUTPUT_FOLDER = 'output_06_07_full'
FASTA_FOLDER = f'{OUTPUT_FOLDER}/fasta'
INTERSECTED_TRANSCRIPTS_FILENAME = f'{OUTPUT_FOLDER}/intersected_transcripts.json'
UTR_STATS_FILENAME = f'{OUTPUT_FOLDER}/utr_statistics.json'
ALL_TRANSCRIPTOME_UTRS_FILENAME = f'{OUTPUT_FOLDER}/all_transcriptome_UTRs.fasta'

NO_HIT_DESC = '***no hit found***'
CATEGORIES = {'NO_HITS': 1, 'FULL': 3, 'PARTIAL': 2}
STOP = 'TAA'


def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Strand lengths are not equal!")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def parse_blast_file(path):
    int_fields = ['qlen', 'alen', 'qstart', 'qend', 'sstart', 'send']
    float_fields = ['pident', 'alen_qlen']
    result = {}
    headers = []
    with open(BLAST_FILENAME) as f:
        for i, line in enumerate(f):
            line = line.strip()

            if line == '':
                continue

            if i == 0:
                headers = line.split('\t')
                continue

            splitted = line.split('\t')
            seq_id = splitted[0].split(' ')[0]

            if seq_id not in result:
                result[seq_id] = None

            if NO_HIT_DESC not in line:
                result[seq_id] = dict([[headers[i], s] for i, s in enumerate(splitted)])

                for header, val in result[seq_id].items():
                    if header in int_fields:
                        result[seq_id][header] = int(val)
                    elif header in float_fields:
                        result[seq_id][header] = float(val)

                    if header in ['qframe', 'sframe']:
                        if val == '1':
                            result[seq_id][header] = '+'
                        elif val == '-1':
                            result[seq_id][header] = '-'
                        else:
                            raise(
                                Exception(f'Blast frame can be only 1/-1 (val {val}), hit {line}'))


    return result


def contains(hit_coords, stop_coords):
    hit_coords = sorted(hit_coords)
    stop_coords = sorted(stop_coords)
    return hit_coords[0] <= stop_coords[0] and hit_coords[1] >= stop_coords[1]


REVERSE_TABLE = {
        ("+", "+", "+"): { "reverse_transcript": False, "reverse_genome": False,
                              "sstart": "sstart", "send": "send", "qstart": "qstart", "qend": "qend"},
        ("-", "+", "+"): { "reverse_transcript": True, "reverse_genome": False,
                              "sstart": "sstart", "send": "send", "qstart": "qend", "qend": "qstart"},
        ("+", "-", "+"): { "reverse_transcript": True, "reverse_genome": False,
                              "sstart": "sstart", "send": "send", "qstart": "qend", "qend": "qstart"},
        ("+", "+", "-"): { "reverse_transcript": True, "reverse_genome": True,
                              "sstart": "send", "send": "sstart", "qstart": "qend", "qend": "qstart"},
        ("-", "+", "-"): { "reverse_transcript": False, "reverse_genome": True,
                              "sstart": "send", "send": "sstart", "qstart": "qstart", "qend": "qend"},
        ("+", "-", "-"): { "reverse_transcript": False, "reverse_genome": True,
                              "sstart": "send", "send": "sstart", "qstart": "qstart", "qend": "qend"},

        ("-", "-", "+"): { "reverse_transcript": False, "reverse_genome": False,
                              "sstart": "sstart", "send": "send", "qstart": "qstart", "qend": "qend"},
        ("-", "-", "-"): { "reverse_transcript": True, "reverse_genome": True,
                              "sstart": "send", "send": "sstart", "qstart": "qend", "qend": "qstart"}
    }

def which_reverse_needed(tran_direction, genome_direction, gene_direction):
    key = (tran_direction == "+", genome_direction == "+", gene_direction == "+")
    return REVERSE_TABLE[key]

def main():
    total_counts = {
        'total_transcript': 0,
        'with_hits_on_genome': 0,
        'PARTIAL':{
            'with_more_than_one_gene': 0,
            'without_genes': 0,
            'with_one_gene': 0,

            'starts_with_STOP': 0,
            'starts_not_with_STOP': 0,

            'where_UTR_EQUAL': 0,
            'where_UTR_NOT_EQUAL': 0
        },
        'without_hits_on_genome': 0
    }
    transcript = {}
    blast = {}

    for seq in SeqIO.parse(TRANSCRIPT_FILENAME, 'fasta'):
        transcript[seq.id] = seq

    blast = parse_blast_file(BLAST_FILENAME)

    stats = {CATEGORIES['NO_HITS']: [],
             CATEGORIES['FULL']: [],
             CATEGORIES['PARTIAL']: []}

    for seq_id, rec in transcript.items():
        if seq_id not in blast:
            raise(Exception(f"Something wrong, {seq_id} not in hit file"))

        total_counts['total_transcript'] += 1

        hit = blast[seq_id]

        # if no hits
        if not hit:
            total_counts['without_hits_on_genome'] += 1
            stats[CATEGORIES['NO_HITS']].append(seq_id)
            continue
        else:
            total_counts['with_hits_on_genome'] += 1

        # if full coverage

        if hit['pident'] == 100.0 and hit['qlen'] == hit['alen']:
            stats[CATEGORIES['FULL']].append(seq_id)
            continue

        stats[CATEGORIES['PARTIAL']].append(seq_id)

    for name, cat_id in CATEGORIES.items():
        print(f"{name}: {len(stats[cat_id])}")
    print()

    # getting statistics of transcripts with partial hits

    genome = SeqIO.to_dict(SeqIO.parse(GENOME_FILENAME, "fasta"))
    annotation_by_contig_id = dict([[e, []] for e in genome.keys()])

    with open(ANNOTATION_FILENAME) as f:
        for line in f:
            line = line.strip()
            if line == '':
                continue

            splitted = line.split('\t')

            record = {}

            record['contig_id'] = splitted[0]
            record['gene_id'] = splitted[-1].split(';')[0].split('=')[-1]
            record['left'] = int(splitted[3])
            record['right'] = int(splitted[4])
            record['direction'] = splitted[6]

            if record['direction'] == '-':
                interval = I.closed(
                    record['left']-3, record['left']-1)
            elif record['direction'] == '+':
                interval = I.closed(
                    record['right']+1, record['right']+3)
            else:
                raise(Exception('Wrong annotation line'))

            record['gene_stop_crds'] = [interval.lower, interval.upper]

            # fucking continue man
            if record['contig_id'] == 'Ctg28_length_81690':
                continue

            annotation_by_contig_id[record['contig_id']].append(record)

    transcripts_with_stops = []
    if exists(INTERSECTED_TRANSCRIPTS_FILENAME):
        with open(INTERSECTED_TRANSCRIPTS_FILENAME) as f:
            transcripts_with_stops = json.load(f)
        total_counts['PARTIAL']['without_genes'] = len(stats[CATEGORIES['PARTIAL']]) - len(transcripts_with_stops)
    else:
        for seq_id in tqdm(stats[CATEGORIES['FULL']]):
            hit = blast[seq_id]
            # hit_interval = I.closed(*sorted([hit['sstart'], hit['send']]))
            hit_coords = [hit['sstart'], hit['send']]

            intersecting_stop_genes = []

            for record in annotation_by_contig_id[hit['sseqid']]:
                stop_interval = record['gene_stop_crds']

                if contains(hit_coords, stop_interval):
                    intersecting_stop_genes.append(record)

            if len(intersecting_stop_genes) > 0:
                data = {}
                data['transcript_id'] = seq_id
                data['hit'] = hit
                data['genes'] = intersecting_stop_genes

                transcripts_with_stops.append(data)
            else:
                total_counts['PARTIAL']['without_genes'] += 1

            with open(INTERSECTED_TRANSCRIPTS_FILENAME, 'w') as out_f:
                json.dump(transcripts_with_stops, out_f, sort_keys=False, indent=4)


    # extracting sequences from transcript & gene for each hit
    utr_stats = {'NOT_EQUAL': {'cnt': 0, 'records': []},
                 'EQUAL': {'cnt': 0, 'records': []}}

    reverse_cases_stats = dict([[e, {'trn_with_STOP': 0, 'trn_without_STOP': 0,
    'gnm_with_STOP': 0, 'gnm_without_STOP': 0 }] for e in REVERSE_TABLE])

    cnt = 0
    for transcript_object in transcripts_with_stops:
        if len(transcript_object['genes']) > 1:
            # TODO: remove it, take the last stop
            total_counts['PARTIAL']['with_more_than_one_gene'] += 1
            continue
        total_counts['PARTIAL']['with_one_gene'] += 1

        stats_record = {}

        main_gene = transcript_object['genes'][0]
        hit = transcript_object['hit']

        gene_direction = main_gene['direction']

        genome_seq = genome[hit['sseqid']].seq
        transcript_seq = transcript[transcript_object['transcript_id']].seq

        # MAGIC here
        rule_key = (hit['qframe'], hit['sframe'], gene_direction)
        rules = REVERSE_TABLE[rule_key]

        sstart = hit[rules['sstart']] - 1
        send = hit[rules['send']] - 1
        qstart = hit[rules['qstart']] - 1
        qend = hit[rules['qend']] - 1

        if rules['reverse_transcript']:
            transcript_seq = transcript_seq.reverse_complement()
            qstart = len(transcript_seq) - qstart - 1
            qend = len(transcript_seq) - qend - 1

        gene_stop = sorted(main_gene['gene_stop_crds'])[0] - 1

        if rules['reverse_genome']:
            genome_seq = genome_seq.reverse_complement()
            gene_stop = len(genome_seq) - gene_stop - 3
            sstart = len(genome_seq) - sstart - 1
            send = len(genome_seq) - send - 1

        transcript_utr_start = qend - (send - gene_stop)
        transcript_utr = transcript_seq[transcript_utr_start:]
        transcript_utr_end = transcript_utr_start + len(transcript_utr)

        genome_utr_start = gene_stop
        genome_utr_end = gene_stop + len(transcript_utr)
        genome_utr = genome_seq[genome_utr_start:genome_utr_end]


        stats_record['genome_seq_id'] = main_gene['contig_id']
        stats_record['genome_utr_start'] = genome_utr_start
        stats_record['genome_utr_end'] = genome_utr_end
        stats_record['transcript_seq_id'] = transcript_object['transcript_id']
        stats_record['transcript_utr_start'] = transcript_utr_start
        stats_record['transcript_utr_end'] = transcript_utr_end
        stats_record['annotated_gene'] = main_gene
        stats_record['transcript_direction'] = hit['qframe']
        stats_record['genome_direction'] = hit['sframe']
        stats_record['gene_direction'] = gene_direction
        stats_record['gen_utr'] = str(genome_utr)
        stats_record['trn_utr'] = str(transcript_utr)

        if len(genome_utr) == len(transcript_utr):
            stats_record['similarity'] = (
                1.0-(hamming_distance(str(genome_utr), str(transcript_utr))/len(genome_utr)))*100
            stats_record['similarity'] = round(
                stats_record['similarity'], 1)
        else:
            stats_record['similarity'] = 0


        if genome_utr.startswith(STOP) and transcript_utr.startswith(STOP):
            total_counts['PARTIAL']['starts_with_STOP'] += 1
        else:
            total_counts['PARTIAL']['starts_not_with_STOP'] += 1



        if genome_utr.startswith(STOP):
            reverse_cases_stats[rule_key]['gnm_with_STOP'] += 1
        else:
            reverse_cases_stats[rule_key]['gnm_without_STOP'] += 1

        if transcript_utr.startswith(STOP):
            reverse_cases_stats[rule_key]['trn_with_STOP'] += 1
        else:
            reverse_cases_stats[rule_key]['trn_without_STOP'] += 1



        if genome_utr == transcript_utr:
            total_counts['PARTIAL']['where_UTR_EQUAL'] += 1
            utr_stats['EQUAL']['records'].append(stats_record)
            utr_stats['EQUAL']['cnt'] += 1
        else:
            total_counts['PARTIAL']['where_UTR_NOT_EQUAL'] += 1
            utr_stats['NOT_EQUAL']['records'].append(stats_record)
            utr_stats['NOT_EQUAL']['cnt'] += 1



    # sorting stats by similarity
    sort_approach = lambda x: x['similarity']

    for key in ['NOT_EQUAL', 'EQUAL']:
        utr_stats[key]['records'] = sorted(
            utr_stats[key]['records'], key=sort_approach)

    # extract data from the code of transcript
    # more_than_one = [e for e in transcripts_with_stops if len(e['genes']) == 5]
    with open(UTR_STATS_FILENAME, 'w') as out_f:
        json.dump(utr_stats, out_f, sort_keys=False, indent=4)


    with open(UTR_STATS_FILENAME+'.csv', 'w') as f:
        delimiter = ','
        header = ["genome_seq_id",
                  "genome_utr_start",
                  "genome_utr_end",
                  "transcript_seq_id",
                  "transcript_utr_start",
                  "transcript_utr_end",
                  "transcript_direction",
                  "genome_direction",
                  "gene_direction",
                  "similarity",
                  "gen_utr",
                  "trn_utr"]
        f.write(f'{delimiter.join(header)}\n')

        for record in utr_stats['NOT_EQUAL']['records']:
            vals = [str(record[h]) for h in header]
            f.write(f'{delimiter.join(vals)}\n')


    transcript_utrs = []

    for record in utr_stats['NOT_EQUAL']['records']:
        similarity_str = str(record["similarity"]).replace(".", "_")
        filename = f'{ similarity_str }_{record["transcript_seq_id"]}'
        with open(f"{FASTA_FOLDER}/{filename}.fasta", 'w') as out_f:
            out_f.write(f">GENOME_UTR\n")
            out_f.write(f"{ record['gen_utr'] }\n")
            out_f.write(f">TRANSCRIPTOME_UTR\n")
            out_f.write(f"{ record['trn_utr'] }\n")

        if record["similarity"] and record["similarity"] < 100:
            rec = SeqRecord(Seq(record['trn_utr']), filename)
            rec.description = ''
            transcript_utrs.append(rec)

    SeqIO.write(transcript_utrs, ALL_TRANSCRIPTOME_UTRS_FILENAME, "fasta")

    pprint(total_counts, sort_dicts=False)
    print()
    pprint(reverse_cases_stats, sort_dicts=False)
    print("Done")


if __name__ == '__main__':
    main()

    # TODO
    # total_transcript: bla
    #   with_hits_on_genome: bla
    #       with_genes: bla
    #           with_STOP_codons: bla
    #               where_UTR_EQUAL: bla
    #               where_UTR_NOT_EQUAL: bla
    #           without_STOP_codons: bla
    #       without_genes: bla
    #   without_hit_on_genome: bla
