# Synopsis
Week two holds most of the preprocessing work. Using ABSD from the Pasteur Institute, and a small sample, this week I learned about the effectiveness
of linClust and RIOT in processing amino acid sequences. 

## Extract ABSD Fasta Sample
From the ABSD database, the first 1000 entries were used as a preliminary sample to process. The fasta file came with standard specifications.
There was the sequence, species, and other pertinent information to the entries.

## Utilize RIOT on Sample
Riot processed the sequences and worked well, especially in terms of efficiency. 

## Utilize LinClust on RIOT Results
LinClust was also able to work very efficiently in terms of clustering the results. The biggest considerations for next week's improvements lie here. 
The nature of the ABSD database is known in terms of the order the entries follow. With no particular sorting to the entries, utilizing the first 1000 entries
resulted in a clustering pattern that was predictable, with 3 clusters highly skewed towards one area or the other. 

## Next Steps
I would recommend next steps to be organizing the fasta entires in ABSD in a meaningful way such that it does provide meaningful results through the RIOT-LinClust pipeline 
built in this week's work. 
