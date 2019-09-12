# Put the original file in input.csv
# Put the output of the taxon finder in input.json
# Run
# It spits the result out in  output.csv



import json
import csv
from dataclasses import dataclass
import argparse

# construct the argument parser and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True,
    help="path and name of the json input resulting from the gnfinder search")
ap.add_argument("-c", "--csv", required=True,
    help="path and name of the corresponding csv to treat")
ap.add_argument("-o", "--output", required=True,
    help="path and name of the treated output ")
args = vars(ap.parse_args())

injson = args["input"]
intab = args["csv"]
outtab = args["output"]


@dataclass
class FoundTaxon:
    """Class for keeping data for each taxon found."""
    start: int
    end: int
    name: str
    classificationPath: str
    dataSourceQuality: str
    dataSourceTitle: str




# Storage for all entries found
indices_and_fields = []


# Grabbing the needed data from the json file
with open(injson, 'r') as j:
    j_parsed = json.loads(j.read())

    for i in j_parsed['names']:

        indices_and_fields.append(FoundTaxon(i['start'], i['end'], i.get('name','NA'), i['verification'].get('classificationPath','NA'),
         i['verification'].get('dataSourceQuality','NA'), i['verification'].get('dataSourceTitle','NA')))

#        print(indices_and_fields)


# Processing the csv and adding the field
beginning = True
last_pos = 0

output_csv = []

# Create output file
with open(outtab, 'w') as csv_file:
    output_csv = csv.writer(csv_file, delimiter='\t', quotechar='"')
    with open(intab, 'r') as csv_file:
        line = csv_file.readline().strip()
        while line:
            row = list(csv.reader([line],delimiter='\t'))[0]
            # Treating the first line to add the new column
            if beginning is True:
                row += ['FoundName']
                row += ['TaxPath']
                row += ['DataSourceTitle']
                beginning = False

            # Storage of all the taxons found for that entry
            fields = []
            paths = []
            datasources = []

            # Look which field has the indice we want and add that to the fields list
            # If the start and end of a field are between the end of the last line and the end of that one
            # we add it to the fields array
            for field in indices_and_fields :
                if (field.start > last_pos) and (field.end < csv_file.tell()):
                    fields += [field.name]
            # Only output calssificationPath if data quality is OK
                #if (field.start > last_pos) and (field.end < csv_file.tell()): #and field.dataSourceQuality == "HasCuratedSources":
                    paths += [field.classificationPath]
                    datasources += [field.dataSourceTitle]


            # Join the taxons with a ;
            row += [';'.join(fields)]
            row += [';'.join(paths)]
            row += [';'.join(datasources)]

            # We now store this last position for the next cycle
            last_pos = csv_file.tell()

            line = csv_file.readline().strip()

            output_csv.writerow(row)

