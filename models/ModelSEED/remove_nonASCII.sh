# This code removes non-ASCII characters from the files downloaded from the ModelSeed
# Ali R. Zomorrodi - Segre Lab @ BU 
# March 24 2015

# Replace <92> with ' (in vi editor do :%s/[\x92]/'/g)
cat ModelSeed_compounds.tsv |sed "s/\x92/'/g" > tmp
mv tmp ModelSeed_compounds.tsv

cat ModelSeed_reactions.tsv |sed "s/\x92/'/g" > tmp
mv tmp ModelSeed_reactions.tsv

# Replace &#39; with '
cat ModelSeed_compounds.tsv |sed "s/&#39;/'/g" > tmp
mv tmp ModelSeed_compounds.tsv

cat ModelSeed_reactions.tsv |sed "s/&#39;/'/g" > tmp
mv tmp ModelSeed_reactions.tsv

# Remove " (in vi editor do :%s/[\xa8]//g
cat ModelSeed_compounds.tsv |sed "s/\xa8//g" > tmp
mv tmp ModelSeed_compounds.tsv

cat ModelSeed_reactions.tsv |sed "s/\xa8//g" > tmp
mv tmp ModelSeed_reactions.tsv

# Remove \xe3 (in vi editor do :%s/[\xe3]//g
cat ModelSeed_compounds.tsv |sed "s/\xe3//g" > tmp
mv tmp ModelSeed_compounds.tsv

cat ModelSeed_reactions.tsv |sed "s/\xe3//g" > tmp
mv tmp ModelSeed_reactions.tsv

# Remove \x80 (in vi editor do :%s/[\x80]//g
cat ModelSeed_compounds.tsv |sed "s/\x80//g" > tmp
mv tmp ModelSeed_compounds.tsv

cat ModelSeed_reactions.tsv |sed "s/\x80//g" > tmp
mv tmp ModelSeed_reactions.tsv
