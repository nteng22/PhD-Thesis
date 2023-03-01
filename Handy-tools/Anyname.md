# Anyname.sh 
a script to remove part of a file name for all files in a directory <br>
e.g. removing PID_XXX from all fasta files <br>
From @raymondkiu (thanks)

all the lines will essentially be:
```
mv old1 new1
mv old2 new2

```
### Step 1: list all the old names in a txt

```
ls > list1.txt
```

### Step 2: substitute "what you want to remove" with nothing

```
sed 's/PID-XXXX//g' list1.txt > list2.txt

        # Syntax is usually: sed 's/regexp/replacement/g'
        # s = substitute, g = global. All matching occurrences in the line will be replaced.

paste list1.txt list2.txt > list3.txt
awk '{print "mv "$1" "$2}' list3.txt > rename.sh
```

### Step 3 run the script

```
bash ./rename.sh
```
