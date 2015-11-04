# 2d-pharmacophore-search

I am interested in attempting some scaffold hopping type experiments

search2d.py takes a initial model computes its 2d pharmacophore using RDKit. This is then search against a smiles library of compound In practice I have been using the Zinc Clean Leads database but a truncated version is uploaded here for trial (Zinc license prevents distribution of large chunks of Zinc)

Ideally I would parallelise the search but at the moment it is a serial search so is on the slow side. The starting molecule is compared with each member of the database using a Tanimoto similarity on the Pharmacophore fingerprint.

python search2d.py > log.txt

Once we have a list of compounds of similar pharmacophore, brics-scaffold-hop.py takes the compounds fragments them using the BRICS algorithm and then from these generates a new library of compounds. We assess these in a negative design fashion by comparing both similarity to the target parent molecule and to a known non-specific inhibitor. 

python brics-scaffold-hop.py > bricsOut.txt

This script also plots a graph of DeltaSimilarity(Positive design - negative design) vs similarity to positive.

<img src="https://github.com/arthuc01/2d-pharmacophore-search/blob/master/brics-scaffold-assess-test200.png?raw=true" />

New generated compounds generated show enriched delta similarity compared to a search against the Zinc DB alone.


