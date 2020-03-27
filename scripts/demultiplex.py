
# The script searches for barcodes in forward and reverse reads
# that are not at the start of the read, but instead it looks for the bcs
# in the first 30 bases (this can be changed by changing the variable
# 'search_until'.

import sys
import gzip
# from Bio import SeqIO
import time


def find_bcs(readpair, sample_data, search_until):

#    print("checking: \n%s\n%s" %(readpair[1],readpair[5]))

    # Try forward orientation, i.e. foward barcode in forwared read and reverse barcode in reverse read
    for sample in sample_data:

        startindex = -1;
        endindex = -1;
        forw = sample_data[sample]['bcs'][0].upper()
        reve = sample_data[sample]['bcs'][1].upper()
#	print "trying: %s\t%s" %(forw,reve)
        if forw in readpair[1][:search_until]:
            startindex = readpair[1].index(forw)
#                print "found forward %s -> %s" %(forw, startindex)
            if reve in readpair[5][:search_until]:
                endindex = readpair[5].index(reve)
#                       print "found reverse %s -> %s" %(reve,endindex)
                break

    if startindex >= 0 and endindex >= 0:
#    	print "forward assigned to sample: %s\n" %sample
        readpair[1] = readpair[1][startindex+len(forw):]
        readpair[3] = readpair[3][startindex+len(forw):]
        readpair[5] = readpair[5][endindex+len(reve):]
        readpair[7] = readpair[7][endindex+len(reve):]
        sample_data[sample]['seqs']['R1'].extend(readpair[:4])
        sample_data[sample]['seqs']['R2'].extend(readpair[4:])

    #    	print "\n%s\n%s" %(readpair[1],readpair[5])
        return

    else:

    # Try reverse orientation, i.e. foward barcode in reverse read and reverse barcode in forward read
    #	print "try reverse\n%s\n%s" %(readpair[1],readpair[5])
        for sample in sample_data:
                startindex = -1;
                endindex = -1;
		#assign barcodes in opposite order
                forw = sample_data[sample]['bcs'][1].upper()
                reve = sample_data[sample]['bcs'][0].upper()
    #        	print "trying: %s\t%s" %(forw,reve)
                if forw in readpair[1][:search_until]:
                        startindex = readpair[1].index(forw)
    #                	print "found forward %s -> %s" %(forw, startindex)
                        if reve in readpair[5][:search_until]:
                                endindex = readpair[5].index(reve)
    #                        	print "found reverse %s -> %s" %(reve,endindex)
                                break

    if startindex >= 0 and endindex >= 0:
#        print "reverse assigned to sample: %s\n" %sample
        readpair[1] = readpair[1][startindex+len(forw):]
        readpair[3] = readpair[3][startindex+len(forw):]
        readpair[5] = readpair[5][endindex+len(reve):]
        readpair[7] = readpair[7][endindex+len(reve):]
        sample_data[sample]['seqs']['R1'].extend(readpair[:4])
        sample_data[sample]['seqs']['R2'].extend(readpair[4:])
#    	print "\n%s\n%s" %(readpair[1],readpair[5])
        return
    else:
    #	print "no proper hit\n"
        return readpair
#    	sample_data['invalid']['R1'].extend(readpair[:4])
#    	sample_data['invalid']['R2'].extend(readpair[4:])

#    print "\n%s\n%s" %(readpair[1],readpair[5])


def touch_files():

    for sample in sample_data:
        fh1 = open(target+'/'+sample+'.R1.fastq', 'w')
        fh2 = open(target+'/'+sample+'.R2.fastq', 'w')
    fh1 = open(target+'/invalid.R1.fastq', 'w')
    fh2 = open(target+'/invalid.R2.fastq', 'w')


def write_out(reads=0):

    for sample in sorted(sample_data):
        if len(sample_data[sample]['seqs']['R1']) > 0:
            #fh1 = open(target+'/'+sample+'.R1.fastq','a')
            fh1 = gzip.open(target+'/'+sample+'.R1.fastq.gz', 'a')
            #fh2 = open(target+'/'+sample+'.R2.fastq','a')
            fh1 = gzip.open(target+'/'+sample+'.R1.fastq.gz', 'a')
    #	    for seq in sample_data[sample]['seqs']:
    #	    	    fh.write(seq+'\n')
            for i in range(len(sample_data[sample]['seqs']['R1'])):
                fh1.write(sample_data[sample]['seqs']['R1'][i]+"\n")
                fh2.write(sample_data[sample]['seqs']['R2'][i]+"\n")
            fh1.close()
            fh2.close()
            sample_data[sample]['count'] += len(sample_data[sample]['seqs']['R1'])
            sample_data[sample]['seqs']['R1'] = []
            sample_data[sample]['seqs']['R2'] = []
            if reads:
                print("%s\t%i read pairs (%.2f %%)" %(sample, sample_data[sample]['count']/4, (float(sample_data[sample]['count'])/4)/reads*100))
#            else:
        #	    print "no valid reads found for sample '%s'" %sample

    if len(invalid_recs['R1']) > 0:
        fh1 = open(target+'/invalid.R1.fastq','a')
        fh2 = open(target+'/invalid.R2.fastq','a')
        for i in range(len(invalid_recs['R1'])):
            fh1.write(invalid_recs['R1'][i]+"\n")
            fh2.write(invalid_recs['R2'][i]+"\n")
        invalid_recs['count'] += len(invalid_recs['R1'])
        invalid_recs['R1'] = []
        invalid_recs['R2'] = []

        fh1.close()
        fh2.close()
        if reads:
            print("\ninvalid\t%i read pairs (%.2f %%)\n" %(invalid_recs['count']/4, (float(invalid_recs['count'])/4)/reads*100))

    if reads:
        print("total number of read pairs processed: %i" %reads)

def process_pairs(f1, f2, sample_data, invalid_recs, search_until):
    """Interleaves two (open) fastq files.
    """
    count = 0
    touch_files()
    while True:
        lines = []
        line = f1.readline()
        if line.strip() == "":
            break
        lines.append(line.strip())

        for i in range(3):
            lines.append(f1.readline().strip())

        for i in range(4):
            lines.append(f2.readline().strip())

        temp = find_bcs(lines, sample_data, search_until)
        if temp:
            invalid_recs['R1'].extend(temp[:4])
            invalid_recs['R2'].extend(temp[4:])
        count+=1
        if (count % 100000) == 1:
            print("["+time.strftime("%c")+"] - %i read pairs processed" %(count/2*2))
            write_out(0)
    return count




search_until = 30

if not len(sys.argv) == 5:
        print("Expecting 4 arguments\n")
        sys.exit()


print(sys.argv[1])
file1 = sys.argv[2]
file2 = sys.argv[3]
target = sys.argv[4]

print(sys.version)

fh = open(sys.argv[1],'r')
sample_data = {}

for l in fh:
    print(l)
    cols = l.strip().split("\t")
    sample = cols[1]
    bcs = cols[2].split(":")
    sample_data[sample] = {'count': 0, 'bcs':[], 'seqs':{ 'R1': [], 'R2': []}}
    sample_data[sample]['bcs'] = bcs


readcount = 0
invalid_recs = {'count':0, 'R1':[], 'R2':[]}

if file1[-2:] == "gz":
    import gzip
    with gzip.open(file1, 'rt') as f1:
        with gzip.open(file2, 'rt') as f2:
            readcount = process_pairs(f1, f2, sample_data, invalid_recs, search_until)
else:
    with open(file1) as f1:
        with open(file2) as f2:
            readcount = process_pairs(f1, f2, sample_data, invalid_recs, search_until)
f1.close()
f2.close()

write_out(readcount)
