import ViennaRNA as RNA

# 1. Toy Sequence (10 nt)
toy_seq = "GCGCAUACGC"
toy_struct, toy_mfe = RNA.fold(toy_seq)

# 2. Official Moderna Sequence (44 nt)
official_seq = "GGAGCAAAACUUGUCGAUUGAGAACAAAAUACAGAAUUUGCUUG"
official_struct, official_mfe = RNA.fold(official_seq)

print("==================================================")
print("         MILESTONE 1: CLASSICAL BENCHMARK         ")
print("==================================================")
print(f"1. TOY SEQUENCE ({len(toy_seq)} nt):")
print(f"   Sequence  : {toy_seq}")
print(f"   Structure : {toy_struct}")
print(f"   Energy    : {toy_mfe:.2f} kcal/mol\n")

print(f"2. MODERNA OFFICIAL SEQUENCE ({len(official_seq)} nt):")
print(f"   Sequence  : {official_seq}")
print(f"   Structure : {official_struct}")
print(f"   Energy    : {official_mfe:.2f} kcal/mol")
print("==================================================")
