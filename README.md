# materialproject_ml
Project on clustering all known ternary compounds (compounds containing 3 different elements) found in nature. This project involves: 
1) Scraping structures from various databases such as ICSD and Materials Project \n
2) Cleaning the data and converting the structure objects to Pickle/Json formats
3) Computing representative fingerprint functions to encode spatial distribution of properties such as atomic number, presence of atoms, susceptibility, oxidation number (this involved a lot of further processing) etc
4) Clustering these fingerprints to find which structures were "Similar" so as to detemine the fundamental kinds of compounds that occur
5) Checking if clusters resemble known compounds types 

Please look at the notebooks and src folders to find the intermediate results and code.
