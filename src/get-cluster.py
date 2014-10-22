#formatting CSV file for WEKA usage
#extracts the cluster from the first column of the TUD data files and appends it onto the end

import csv

#writes to a new file, to preserve old data
#filepath = '/Users/edwardwilliams/Documents/tango/data/kmer_counts/Hatfull_60_subset/Hatful_60_dev_k4.csv'
#newfile = '/Users/edwardwilliams/Documents/tango/data/kmer_counts/Hatfull_60_subset/Hatful_60_dev_k4_with_clusters.csv'
filepath = '/Users/edwardwilliams/Documents/tango/data/kmer_counts/663_phage/all_dev_k4.csv'
newfile = '/Users/edwardwilliams/Documents/tango/data/kmer_counts/663_phage/all_dev_k4_with_clusters.csv'

with open(filepath) as csvfile:
	with open(newfile, 'w') as new:
		writer = csv.writer(new)
		reader = csv.reader(csvfile)
		for row in reader:
			if row[0] != "ID": #i.e. it's not the first row
				cluster = row[0].split("(")[1]#pulls out whatever is between the first pair of parenthesis
				finalcluster = cluster.split(")")[0]
				print finalcluster
				row.append(finalcluster)
				writer.writerow(row)
			else:
				row.append("Cluster")
				writer.writerow(row) 


