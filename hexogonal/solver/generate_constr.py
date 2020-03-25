def generate_constr(filename, length):
    constr_file = open(filename , 'w')
    constr_file.write("12" + '\n')
    constr_file.write("0 xy 0.0 0.0" + '\n')
    constr_file.write("13 xy 0.0 6.20546" + '\n')
    constr_file.write("26 xy 0.0 7.997262150038588" + '\n')
    constr_file.write("39 xy 0.0 14.20272215003859" + '\n')
    constr_file.write("52 xy 0.0 15.994524300077176" + '\n')
    constr_file.write("65 xy 0.0 22.199984300077176" + '\n')

    constr_file.write("12 xy " + length + " 0.0" + '\n')
    constr_file.write("25 xy " + length + " 6.20546" + '\n')
    constr_file.write("38 xy " + length + " 7.997262150038588" + '\n')
    constr_file.write("51 xy " + length + " 14.20272215003859" + '\n')
    constr_file.write("64 xy " + length + " 15.994524300077176" + '\n')
    constr_file.write("77 xy " + length + " 22.199984300077176" + '\n')
    constr_file.close()



if __name__ == "__main__":
    length = "35.82426"
    number = 21
    filename = "constr3x6_" + str(number) + ".txt"
    generate_constr(length, filename)