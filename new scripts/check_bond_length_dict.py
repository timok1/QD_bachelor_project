import bonding_distances as bond_dis
import helper_functions as hf


def check_in_dict(el):
    check_bool = False
    for entry in bond_dis:
        if entry == el:
            print(el + " is in dictionary")
            return check_bool
    return check_bool


def enter_dis2dict(el):
    for entry in bond_dis:
        val = float(input("Add bonding length between " + el + " and " + entry + ": "))


if __name__ == "__main__":
    check = True
    while check:
        el = input("Enter element to be checked, or type 'all': ")
        if el == 'all':
            print(dir(bond_dis))
            for entry in bond_dis.builtins:
                print(entry)
        else:
            check = check_in_dict(el)
    enter_val = hf.y2true(input(el + " is not currently in dictionary. Do you want to add the bonding distances? (y/n): "))
    if enter_val:
        enter_dis2dict(el)
