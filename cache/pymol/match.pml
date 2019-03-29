load E:/CSRC/Nepre/NEPRE-DESKTOP/cache/5fjl.pdb, na
load E:/CSRC/Nepre/NEPRE-DESKTOP/cache/5fjl.pdb, dcy
super na, dcy, object=aln
python
pairs = cmd.get_raw_alignment('aln')
print pairs[0], pairs[-1]
print pairs[0][1][1], pairs[-1][1][1], pairs[0][1][1]-pairs[-1][1][1]
cmd.select("sel1","na and index %d-%d"%(pairs[0][0][1],pairs[-1][0][1]))
cmd.save("../cache/files/match_segment.pdb","sel1")
python end