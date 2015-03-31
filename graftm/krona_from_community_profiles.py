
import csv
import subprocess
import tempfile

class OtuTable:
    def __init__(self, sample_name):
        self.sample_name = sample_name
        self.sampleCounts = {} #taxonomy to count info


class KronaBuilder:
    def otuTablePathListToKrona(self, otuTablePaths, outputName):
        otuTables = []
        for table in self.parseOtuTable(otuTablePaths):
            otuTables.append(table)

        self.runKrona(otuTables, outputName)

    def parseOtuTable(self, otuTablePath):
        data = csv.reader(open(otuTablePath), delimiter="\t")

        # Parse headers (sample names)
        fields = data.next()
        if len(fields) < 3: raise "Badly formed OTU table %s" % otuTablePath
        tables = []
        for i in range(len(fields)-2):
            table = OtuTable(fields[i+1])
            tables.append(table)

        # Parse the data in
        taxonomyColumn = len(fields)-1
        for row in data:
            for i in range(len(fields)-2):
                taxonomy = row[taxonomyColumn]
                tables[i].sampleCounts[taxonomy] = row[i+1]

        return tables

    def runKrona(self, otuTables, outputName):

        # write out the tables to files
        tempfiles = []
        tempfile_paths = []
        for table in otuTables:
            tmps = tempfile.mkstemp('','CommunityMkrona')
            tmp = tmps[1]
            out = open(tmp,'w')
            tempfiles.append(out)
            tempfile_paths.append(tmp)
            for taxonomy, count in table.sampleCounts.iteritems():
                tax = "\t".join(taxonomy.split(';'))
                out.write("%s\t%s\n" % (count,tax))
            out.close()

        cmd = ["ktImportText",'-o',outputName]
        for i, tmp in enumerate(tempfile_paths):
            cmd.append(','.join([tmp,otuTables[i].sample_name]))

        # run the actual krona
        subprocess.check_call(' '.join(cmd) + ' 1>/dev/null ', shell=True)

        # close tempfiles
        for t in tempfiles:
            t.close()



