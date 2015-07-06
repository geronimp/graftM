

class HMMreader:
    def __init__(self, table):
        self.entries = {}
        self.type    = [x.split()[2] for x in open(table, 'r') if x.startswith('# Program:')][0]
        table        = [x.rstrip().split() for x in open(table, 'r') if not x.startswith('#')]
        for entry in table:
            if entry[0] in self.entries:
                if self.type == "hmmsearch":    

                    if float(self.entries[entry[0]][15]) > float(entry[15]):
                        self.entries[entry[0]][15]=entry[15]
                    if float(self.entries[entry[0]][16]) < float(entry[16]):
                        self.entries[entry[0]][16]=entry[16]
                    
                elif self.type == "nhmmer":
                    if float(self.entries[entry[0]][13]) > float(entry[13]):
                        self.entries[entry[0]]=entry
            else:
                self.entries[entry[0]]=entry
        
    def names(self):
        return self.entries.keys()

    def evalue(self, entry):
        if self.type == "nhmmer":
            return float(self.entries[entry][12])
        elif self.type == "hmmsearch":
            return float(self.entries[entry][6])
        else: return None

    def name(self, entry):
        return self.entries[entry][0]

    def hmm_len(self, entry):
        if self.type == "nhmmer":
            return str(self.entries[entry][10])
        elif self.type == "hmmsearch":
            return str(self.entries[entry][5])
        else: return None

    def seq_len(self, entry):
        if self.type == "nhmmer":
            return str(self.entries[entry][10])
        elif self.type == "hmmsearch":
            return str(self.entries[entry][2])
        else: return None

    def hmmfrom(self, entry):
        if self.type == "nhmmer":
            return str(self.entries[entry][4])
        elif self.type == "hmmsearch":
            return str(self.entries[entry][15])
        else: return None

    def hmmto(self, entry):
        if self.type == "nhmmer":
            return str(self.entries[entry][5])
        elif self.type == "hmmsearch":
            return str(self.entries[entry][16])
        else: return None

    def aln_len(self, entry):
        if self.type == "nhmmer":
            tofrom=[float(self.entries[entry][5]),float(self.entries[entry][4])]
            return max(tofrom)-min(tofrom)
        elif self.type == "hmmsearch":
            tofrom=[float(self.entries[entry][15]),float(self.entries[entry][16])]
            return max(tofrom)-min(tofrom)
        else: return None

    def alifrom(self, entry):
        if self.type == "nhmmer":
            return str(self.entries[entry][6])
        elif self.type == "hmmsearch":
            return str(self.entries[entry][17])
        else: return None

    def alito(self, entry):
        if self.type == "nhmmer":
            return str(self.entries[entry][7])
        elif self.type == "hmmsearch":
            return str(self.entries[entry][18])
        else: return None

    def strand(self, entry):
        if self.type == "nhmmer":
            return self.entries[entry][11]
        elif self.type == "hmmsearch":
            return None
        else: return None

    def bit(self, entry):
        if self.type == "nhmmer":
            return float(self.entries[entry][13])
        elif self.type == "hmmsearch":
            return float(self.entries[entry][7])
        else: return None
