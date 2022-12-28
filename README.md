
# Git notes:

Always do the following, when first starting out:
0. git checkout master
1. git pull
2. git checkout local
3. git merge master

4. ---> Do your work now, when done, go to step 5 <---

5. git add (files you want to add) or . to add all files
6. git commit -m "message here"
7. git checkout master
8. git merge local
9. git push 

If you want to cancel local changes, but save them for possible re-use later:
1. git stash
This command will reset changes to all files, but also saves them in case you need them later. 

If you need the files later you can always do:
2. git stash pop 
