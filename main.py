# import gc, Datastructure
# import sys

# #boolean functions

# # dna = input('Enter your DNA Squance, please:')
# # frame = input('framing number:')
# # CodonFound =''
# # def hasStopCodon(dna,frame = 0):
# #   CodonFound = False
# #   StopCodons = ['tag', 'tag','taa']
# #   for i in range(frame,len(dna),3):
# #     codon=dna [i:i+3].lower()
# #     if codon in StopCodons :
# #       CodonFound = True
# #       break

# #   return CodonFound

# # print (hasStopCodon ('tagggactagagtacccaaca'))

# # def rev (seq):
# #   seq = reverse(seq)
# #   seq = complement(seq)
# #   return seq

# # print(rev ('ccaaggaaggaa'))

# # dna = 'agtagtaa'
# # dna[0:2:3]
# # dna [::-1]

# # print (dna [0:2:3] , dna [::-1])
# # def complement (dna):
# #   """ return the complementary sequence string."""
# #   basecomplement = { 'a':'t', 't':'a', 'g':'c', 'c':'g', 'A':'T','C':'G','G':'C','T':'A','N':'N'}
# #   letters = list(dna)
# #   letters = [basecomplement[base] for base in letters]
# #   return ''.join(letters)

# # print (complement ('agtagtaa'))
# # import random
# # def create_dna(n, alphabet='acgt'):
# #     return ''.join([random.choice(alphabet) for i in range(n)])

# # dna = create_dna(10000)
# # print(dna)

# # def count1(dna, base):
# #   i = 0
# #   for c in dna:
# #       if c == base:
# #     i += 1
# #   return i

# # def count2(dna, base):
# #   i = 0
# #   for j in range(len(dna)):
# #       if dna[j] == base:
# #     i += 1
# #   return i

# # def count3(dna, base):
# #   match = [c == base for c in dna]
# #   return sum(match)

# # def count4(dna, base):
# #   return dna.count(base)

# # def count5(dna, base):
# #   return len([i for i in range(len(dna)) if dna[i] == base])

# # def count6(dna,base):
# #   return sum(c == base for c in dna)

# # f=open('','r')
# # f=open('', 'w')
# # f=open('', 'a')

# # try:
# #   f=open('myfile')
# # except FileNotFoundError:
# #   print('File not found')

# # f.read('myfile')
# # f.seek(0)
# # for line in f:
# #   print(line)
# # f.readline()
# # f.write('string')

# # f.close()

# # try:
# #   f= open(myfile.fa)
# # except FileNotFoundError:
# #   print('File not found')

# # seqs = {}
# # for line in f:
# #   line=line.rstrip()
# #   if line [0] ==">":
# #     words=line.split()
# #     name=words [0] [1:]
# #     seqs[name]=''
# #   else:
# #     seqs [name] = seqs [name] + line

# # print(seqs)

# # close (myfile.fa)

# # def compute(n,x,y) :
# #   if n==0 : return x
# #   return compute(n-1,x+y,y)

# # print (compute(2,4,6))

# # dna ="TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG "

# # from Bio.Blast import NCBIWWW
# # from Bio.Blast import NCBIXML

# # sequence_data = dna

# # result_handle = NCBIWWW.qblast("blastn", "nt", sequence_data)

# # blast_records = NCBIXML.parse(result_handle)
# # blast_record = next(blast_records)  # Assuming you want the first hit

# # for alignment in blast_record.alignments:
# #     for hsp in alignment.hsps:
# #         print("Sequence:", alignment.title)
# #         print("Length:", alignment.length)
# #         print("E-Value:", hsp.expect)
# #         # Additional details like score, identities, etc.

# # from Bio.Seq import Seq
# # from Bio.Seq import translate

# # dna_sequence = "TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAGN"

# # protein_sequence = translate(dna_sequence)

# # print(protein_sequence)

fasta_file = "dna2.fasta"


def count_fasta_records(file_path):
  """Counts the number of records in a FASTA file."""

  record_count = 0
  with open(file_path, 'r') as file:
    for line in file:
      if line.startswith(">"):
        record_count += 1

  return record_count


num_records = count_fasta_records(fasta_file)
print("Number of records:", num_records)


def analyze_fasta_lengths(file_path):
  """Analyzes sequence lengths in a FASTA file."""

  sequences = {}  # Store sequence identifiers and lengths
  current_id = None

  with open(file_path, 'r') as file:
    for line in file:
      if line.startswith(">"):
        current_id = line[1:].split()[0]  # Extract identifier
        sequences[current_id] = 0
      else:
        sequences[current_id] += len(line.strip())

  longest_seq = max(sequences, key=sequences.get)
  shortest_seq = min(sequences, key=sequences.get)

  return sequences, longest_seq, shortest_seq


all_lengths, longest_id, shortest_id = analyze_fasta_lengths(fasta_file)

# print("Sequence lengths:", all_lengths)
# print("Longest sequence ID:", longest_id, "Length:", all_lengths[longest_id])
# print("Shortest sequence ID:", shortest_id, "Length:", all_lengths[shortest_id])

from Bio import SeqIO


def get_reading_frame(dna_sequence, frame):
  """Returns the specified reading frame of a DNA sequence."""
  start = frame - 1
  return dna_sequence[start::3]


def find_orfs(reading_frame_sequence):
  """Finds all ORFs within a reading frame."""
  start_codon = "ATG"
  stop_codons = ["TAA", "TAG", "TGA"]
  orfs = []
  start_index = None

  for i in range(0, len(reading_frame_sequence), 3):
    codon = reading_frame_sequence[i:i + 3]
    if codon == start_codon:
      start_index = i
    if codon in stop_codons:
      if start_index is not None:
        orfs.append((start_index + 1, i + 3))  # Positions are 1-indexed
      start_index = None
  return orfs


# File handling

longest_orf_length = 0
longest_orf_seq_id = None
longest_orf_start_frame3 = -1

for record in SeqIO.parse(fasta_file, "fasta"):
  seq_id = record.id
  dna_sequence = str(record.seq)

  # Longest ORF per sequence
  longest_orf_in_seq = 0

  for frame in range(1, 4):
    reading_frame_sequence = get_reading_frame(dna_sequence, frame)
    orfs = find_orfs(reading_frame_sequence)

    for start, end in orfs:
      orf_length = end - start

      # Overall longest ORF
      if orf_length > longest_orf_length:
        longest_orf_length = orf_length
        longest_orf_seq_id = seq_id

      # Longest ORF in frame 3
      if frame == 3 and orf_length > longest_orf_start_frame3:
        longest_orf_start_frame3 = start

      # Longest in the specific sequence
      longest_orf_in_seq = max(longest_orf_in_seq, orf_length)

  # Output for the specific sequence
  if seq_id == "gi|142022655|gb|EQ086233.1|16":
    print(f"Length of longest ORF in sequence {seq_id}:", longest_orf_in_seq)

# Overall Outputs
# print("Starting position of longest ORF in reading frame 3:", longest_orf_start_frame3)
# print("Length of the longest ORF (any sequence, any frame):", longest_orf_length)
# print("Sequence identifier containing the longest ORF:", longest_orf_seq_id)


def find_repeats(sequences, n):
  """Finds all repeats of length n in a list of sequences."""
  repeats = {}
  for seq in sequences:
    for i in range(len(seq) - n + 1):
      repeat = seq[i:i + n]
      repeats[repeat] = repeats.get(repeat, 0) + 1  # Count occurrences
  return repeats


# File Handling
sequences = [str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]

# Specify the repeat length
repeat_length = 3

# Find repeats
repeats = find_repeats(sequences, repeat_length)

# Find the most frequent repeat
most_frequent_repeat = None
highest_count = 0

for repeat, count in repeats.items():
  if count > highest_count:
    most_frequent_repeat = repeat
    highest_count = count

# Output
# print(f"Repeats of length {repeat_length} and their occurrences:")
for repeat, count in repeats.items():
  print(f"{repeat}: {count}")

# print(f"\nMost frequent repeat of length {repeat_length}:")
# print(f"{most_frequent_repeat}: {highest_count}")

all_orfs = []  # Store details of all ORFs

for record in SeqIO.parse(fasta_file, "fasta"):
  seq_id = record.id
  dna_sequence = str(record.seq)

  # Longest ORF per sequence and track all ORFs for subsequent analysis
  longest_orf_in_seq = 0

  for frame in range(1, 4):
    reading_frame_sequence = get_reading_frame(dna_sequence, frame)
    orfs = find_orfs(reading_frame_sequence)

    for start, end in orfs:
      orf_length = end - start
      all_orfs.append((seq_id, frame, start, orf_length))  # Track all ORFs

      # ... (Existing longest ORF tracking logic remains the same)

# ... (Functions: find_repeats)

# ---Answers to Your Questions---
# 1. Starting position of the longest ORF in reading frame 3
longest_orf_frame3 = max((orf for orf in all_orfs if orf[1] == 3),
                         key=lambda x: x[3])
print(
    f"Starting position of longest ORF in reading frame 3: {longest_orf_frame3[2]}"
)

# 2. Length of longest ORF (any sequence, any forward reading frame)
longest_orf_any_forward = max((orf for orf in all_orfs if orf[1] in [1, 2, 3]),
                              key=lambda x: x[3])
print(
    f"Length of the longest ORF (any sequence, any forward frame): {longest_orf_any_forward[3]}"
)

# 3. Longest ORF in the specific sequence
target_id = "gi|142022655|gb|EQ086233.1|16"
longest_in_target = max((orf for orf in all_orfs if orf[0] == target_id),
                        key=lambda x: x[3])
print(f"Length of longest ORF in sequence {target_id}: {longest_in_target[3]}")

# 4. Most frequent repeat of length 6
repeat_length = 6
repeats = find_repeats(sequences, repeat_length)
most_frequent_6 = max(repeats, key=repeats.get)
print(
    f"Most frequent repeat (length 6): {most_frequent_6} (occurs {repeats[most_frequent_6]} times)"
)

# 5. Repeats of length 12 and their frequency
repeat_length = 12
repeats_12 = find_repeats(sequences, repeat_length)
max_count = max(repeats_12.values())
most_frequent_12 = [
    repeat for repeat, count in repeats_12.items() if count == max_count
]
print(f"Repeats of length 12 with maximum occurrences:")
for repeat in most_frequent_12:
  print(f"{repeat} (occurs {max_count} times)")

# 6. Most frequent from a given set of repeats
target_repeats = ["CGCGCCG", "CATCGCC", "TGCGCGC"]
all_repeats = find_repeats(sequences, 7)  # Adjust length to 7
max_count = 0
most_frequent = None
for repeat in target_repeats:
  count = all_repeats.get(repeat, 0)  # Get count with default 0 if not found
  if count > max_count:
    max_count = count
    most_frequent = repeat
print(
    f"Most frequent repeat from the list: {most_frequent} (occurs {max_count} times)"
)