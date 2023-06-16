file=open("perform_roc.txt", "r")
lst=[]
fpr_list1=[]
fpr_list=[]
tpr_list=[]
tpr_list1=[]
for line in file:
    lst.append(line.strip().split("["))
for i in range(len(lst)):
    fpr_list1=lst[0][1].split(",")
    tpr_list1=lst[1][1].split(",")
fpr_list1.pop(-1)
tpr_list1.pop(-1)
for el in fpr_list1:
    fpr_list.append(float(el))
for element in tpr_list1:
    tpr_list.append(float(element))
print(fpr_list)
print(tpr_list)

import matplotlib.pyplot as plt
plt.xlabel("False Positive Rate", fontsize=14)
plt.ylabel("True Positive Rate", fontsize=14)
plt.title("Receiver Operating Characteristic Curve", fontsize=14)
plt.plot(fpr_list, tpr_list, color='purple', linewidth=2)