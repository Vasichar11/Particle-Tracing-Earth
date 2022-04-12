import csv 


negatives=[line.strip() for line in open('negative_aeqs.txt')]

for i in negatives:
    with open('dstr_before_and_after_WPI.csv', "r") as f1:
        writer = csv.writer(f1)
        for line in f1:
            if( line.split(",")[0] == i ):  
                #writer.writerow(line.split(","))
                print(line.split(","))




