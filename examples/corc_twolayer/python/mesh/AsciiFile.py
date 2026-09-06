import os
import re

class AsciiFile:
# --------------------------------------------------------------------------

    def __init__(self, Path = "" ):

        ## The path that is associated with this file
        self.Path = Path

        ## The buffer contains all lines in this ascii file
        self.Buffer = []

        # try to load data from file
        if self.Path != "":

            # Test if file Exists
            if not os.path.exists(Path):
                raise Exception('File ' + str(Path) + ' does not exist')

            # load the file from the buffer
            File = open(Path, 'r',encoding='utf-8')

            # add line to file
            for Line in File:
                self.Buffer.append( Line )

            # close the file instance
            File.close()

# --------------------------------------------------------------------------

    def save(self, Path = '' ):

        if( Path == '' ):
            Path = self.Path

        # save file into output
        File = open(Path, 'w',encoding='utf-8')

        # dump buffer into file
        for Line in self.Buffer:
            File.write(Line + '\n')

        # close file
        print('saving ' + Path)
        File.close()

# --------------------------------------------------------------------------

    def __del__(self):
        # delete the buffer
        self.Buffer.clear()

# --------------------------------------------------------------------------

    def clean_buffer(self):

        n = len(self.Buffer)

        # loop over all lines in buffer
        for k in range(n):
            # copy line
            Line = re.sub('\t', ' ', self.Buffer[k])

            # Tidy up line
            self.Buffer[k] = re.sub(' +', ' ', Line.strip())
