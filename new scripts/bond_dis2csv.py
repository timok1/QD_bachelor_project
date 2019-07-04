import helper_functions as hf
import csv
import pandas as pd


def enter_dis2dict(el):
    length_dict = {}
    with open('bonding_distances.csv', 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            x = input("Enter distance between" + el + " and " + row["Element Name"])
            length_dict[row["Element Name"]] = x
        length_dict[el] = input("Enter distance between" + el + " and " + el)
        csvfile.close()
    csv_input = pd.read_csv('bonding_distances')
    csv_input['']



if __name__ == "__main__":
    check = True
    enter_val = False
    while check:
        with open('bonding_distances.csv') as csvfile:
            reader = csv.DictReader(csvfile)
            el = input("Enter element to be checked, or type 'all' to see all supported elements: ")
            if el == 'all':
                for row in reader:
                    print(row["Element Name"])
            else:
                for row in reader:
                    if row["Element Name"] == el:
                        print(el + " is currently supported")
                    else:
                        enter_val = hf.y2true(input(el + " is not currently in dictionary. Do you want to add the bonding distances? (y/n): "))
                        check = False
        if enter_val:
            enter_dis2dict(el)
    csvfile.close()
