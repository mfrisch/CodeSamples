#!/usr/bin/env python

# Takes as input a project file from the
# European Nucleotide Archive SRA
# and downloads the associated project fastq files

import urllib2
import sys

f = open(sys.argv[1],"r")
table = f.read()
table = table.split('\n')
for i in range(0,len(table)):
	table[i] = table[i].split('\t')
f.close()

# Extract fastq file addresses (Fastq files (ftp))
# Can be merged with for loop above, but O(2n) = O(n) so whatever
ftps = ()
for i in range(1,len(table)-1):
	if ";" in table[i][10]:
		temp = table[i][10].split(';')
		ftps += (temp[0],)
		ftps += (temp[1],)
	else:
		ftps += (table[i][10],)

# Download fastq files
# Courtesy of PabloG and Bjorn Pollex
# at http://stackoverflow.com/questions/22676/
for i in range(0,len(ftps)):
	url = "ftp://"+ftps[i]

	file_name = url.split('/')[-1]
	u = urllib2.urlopen(url)
	f = open(file_name, 'wb')
	meta = u.info()
	file_size = int(meta.getheaders("Content-Length")[0])
	print "Downloading: %s Bytes: %s" % (file_name, file_size)

	file_size_dl = 0
	block_sz = 8192
	while True:
	    buffer = u.read(block_sz)
	    if not buffer:
	        break

	    file_size_dl += len(buffer)
	    f.write(buffer)
	    status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
	    status = status + chr(8)*(len(status)+1)
	    print status,

f.close()