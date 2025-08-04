import ROOT
import sys

filename = sys.argv[1]
treename = sys.argv[2]
branchname = sys.argv[3]

f = ROOT.TFile.Open(filename)
tree = f.Get(treename)

total = 0.0
nentries = tree.GetEntries()
print_interval = max(1, nentries // 10)  # Every 10%

for i, entry in enumerate(tree):
    total += getattr(entry, branchname)

    if i % print_interval == 0 or i == nentries - 1:
        percent = int(100 * i / nentries)
        print(f"Processed {i}/{nentries} entries ({percent}%)")

print(f"\nTotal sum of '{branchname}' = {total}")

