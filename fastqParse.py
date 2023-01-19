# !/usr/bin/env python3
# Name: Your full name
# Group Members: List full names or “None”
'''
@EAS139:136:FC706VJ:2:2104:15343:197393 then the program will output:
Instrument = EAS139
Run ID = 136
Flow Cell ID = FC706VJ
Flow Cell Lane = 2
Tile Number = 2104
X-coord = 15343
Y-coord = 197393 ### Hints:
'''


class FastqString(str):
    ''' Class docstring goes here. Data is stored into an object within the bag'''
    def __init__(self, FASTQ):
        self.FASTQ = FASTQ

    def parse(self):
        ''' Here the inputted FASTQ input is run through every semicolon. A +1 is used to skip the initial
        semicolon to search for the next semicolon. Each segment is printed individually to it's belonging
        section statement. +1 is also used to ignore semicolon when printed as well.
        '''
        instrument = self.FASTQ.find(":")
        run_id = self.FASTQ.find(":", (instrument + 1))
        flow_cell_id = self.FASTQ.find(":", (run_id + 1))
        flow_cell_lane = self.FASTQ.find(":", (flow_cell_id + 1))
        tile_number = self.FASTQ.find(":", (flow_cell_lane + 1))
        x_coord = self.FASTQ.find(":", (tile_number + 1))
        y_coord = self.FASTQ.find(":",(x_coord + 1))
        print("Instrument = " + self.FASTQ[1:instrument])
        print("Run ID = " + self.FASTQ[instrument + 1:run_id])
        print("Flow Cell ID = " + self.FASTQ[run_id + 1:flow_cell_id])
        print("Flow Cell Lane = " + self.FASTQ[flow_cell_id + 1:flow_cell_lane])
        print("Tile Number = " + self.FASTQ[flow_cell_lane + 1:tile_number])
        print("X-coord = " + self.FASTQ[tile_number + 1:x_coord])
        print("Y-coord = " + self.FASTQ[x_coord + 1:-1])
        pass


def main():
    ''' Function docstring goes here.'''
    while (True):
        data = input('Insert FASTQ file data:')
        thisFASTQ = FastqString(data)
        pureData = thisFASTQ.parse()
    pass


main()