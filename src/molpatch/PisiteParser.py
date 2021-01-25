
class PisiteParser:
    """
    Parse pisite files to dict

    ...

    Attributes
    ----------
    pisite_file : str
        file with the pisite data

    """
    def __init__(self, file):
        self.file = file
        self.pisite_dict = self.parse()

    def parse(self):
        """parse a pisite file to a dict

        return
        ------
        dict
            {pisite id: list}
        """
        #open file for reading
        print("opening file ", self.file)
        try:
            infile = open(self.file, 'r')
            start_reading = False
            pisite_dict = {}
            for line in infile.readlines():
                line =  line.rstrip()
                if(start_reading):
                    data = line.strip().split()
                    if len(data) < 4:
                        continue
                    pisite_dict[int(data[0])] = data[1:4]
                elif line.startswith('#residue_no'):
                    start_reading = True

            infile.close()
            return pisite_dict
        except IOError:
            print("Error: Cannot open Pisite file " + self.file + ".")
            return

    def interaction_sites(self):
        """Get total interaction sites

        return
        ------
        int
            total interaction sites
        """
        return sum(1 for x in self.pisite_dict if int(self.pisite_dict[x][2]) == 1)

    def get_data(self):
        """Return interaction dict

        return
        ------
        dict
            {pisite id: list}
        """
        return self.pisite_dict
