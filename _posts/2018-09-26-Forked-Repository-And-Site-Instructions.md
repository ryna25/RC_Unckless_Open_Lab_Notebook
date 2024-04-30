---
layout: post
title: Forking the Putnam Open Lab Notebook for Your GitHub
category: [ Make Your Online Lab Notebook ]
---

### Walk-though of Forking and Modifying this Notebook Site
#### All steps that can be done through the terminal instead of a GUI will be written for doing in the terminal because it is good practice to use the shell

## Forking a Blog Style Github for Notebook Purposes


**Important Links:**  
[Original Jekyll Blog](https://github.com/barryclark/jekyll-now)  
[Maggie's Notebook Website](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/)  
[Maggie's Notebook Github Repository](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook)  
**The Jekyll Now blog/website is the best documented and easiest to modify into a lab notebook. This is why this style notebook is used for making your first Github-based online lab notebook. However, that's not the only blog/notebook you can use! Additionally, if you want to make more than one online lab notebook, you cannot fork the same one twice, and will have to chose a different one (this is the reason Maggie has 3 different ones).**  
[More jekyll themes](https://github.com/topics/jekyll-themes)  
[Another blog style notebook Maggie uses for the Unckless Lab](https://meschedl.github.io/Unckless-Lab-Notebook-Maggie/)


### Steps Specifically for the Jekyll Now Style Notebook from Maggie's Putnam Lab Notebook

1. Open your Terminal application. If you have a mac the name of the app is terminal. If you have a pc there should be an option to open the Command Prompt, or you may have to install putty.
2. `cd Desktop/folder` Get yourself into the folder you want, you can make folders with your finder or in the terminal as in the next step. If your folder is in a different place, like Documents put that in the path.
3. `mkdir Notebook-Directory` Make a directory (same thing as a folder) where you want to put the repository for the Putnam Lab Notebook. This is where all the files will go.
4. `cd Notebook-Directory` Go into the directory you just made.
5. `git init` Now make that directory an empty git repository. Git is a way to do version control, which tracks changes in your files and can help you if previous versions of you.
  - You may get an error here that says invalid active developer path, use this code `xcode - select --install` in your terminal to fix it, you should then be able to use git init.
6. Login to [GitHub](https://github.com/) and go to [Maggie's Putnam Lab Notebook Repository](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook). If you don't have an account, make one. It's free!
7. Click the button in the top right corner that says Fork. It will take everything from this repository and create a copy in your account. Once that is done go into the settings of that repository and change the name. For example mine is MESPutnam_Open_Lab_Notebook. Remember to not use spaces in the name, and it is best to make it descriptive so other's can find it on your Github account.
8. Note: whenever you want to get back to the home of your repository just click the <>Code tab. In that tab click the green button on the right that says Clone or Download. If will give you a link, copy that link.
9. You want to make a clone of this repository on your personal computer and have it connect to GitHub. This way you can change things on your computer, then send them to GitHub as a way of version control. In your command prompt or terminal you should still be in the directory you created for this notebook. You can always use the command `pwd` to check your working directory.
10. `git clone https://github.com/GitHub-user-name/Your-Notebook-Repo-Name.git` The terminal will print out some things saying that it's downloading the repository. Now you have a remote copy of your repo on your desktop, and you should be able to see it in your finder window too.
11. Back in GitHub online, go to the settings tab again. Click on the GitHub Pages section; this is where you tell the site that you want your notebook to be made into a website. Make sure the brand deployment is set to the master branch.
12. A box will come up that says your site is ready to be published/is published and then gives you a link that will look something like this: https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/ but with your username and repo name. You might have to refresh for it to come up. Copy this url.
13. Go back to the <>Code section.  At the top right it says "Maggie's Putnam Lab Open Lab Notebook" and the link to Maggie's notebook. Click edit, rename the description to whatever you want, and paste the link you got from the settings page. This makes it so that people can click the link directly when they get to your repository and it will take them to the actual website of your notebook (which at this point is still a copy of Maggie's, but you will change that later).
14. On your computer finder go to the folder where your repository is. Open the file called _config.yml in your preferred text editor. Maggie uses [Atom](https://atom.io/) but that will soon be depreciated. Other text editors are [Visual Studio Code](https://code.visualstudio.com/), [Notepad++](https://notepad-plus-plus.org/), [Sublime Text](https://www.sublimetext.com/), and many others that you can look up. There are many that are good free, so there is no need to purchase one.
15. The _config.yml file is basically a map for your site. Change all of these to personalize the site:
    - Name
    - Description
    - Avatar: This puts in an image you want, if you want just the profile picture on your GitHub account already you can right click on it to and copy the image link. If you want another one you can specify the path to that image.
    - All the social media: This is then put at the footer of all your pages. If you don't have a certain social site then leave it blank
16. You also need to change the lines 44 and 50 in the _config.yml file. Change url: on 44 to equal "http://your-username.github.io" and delete what is after the words baseurl: on line 50
16. It shouldn't be necessary to change any of the other sections, if you mess those up your site won't publish. Save this file. However, that only saves it to your computer. We have to add it to Git and send it to GitHub before it will change anything online. This is part of version control where changes are tracked and your repository will live both on your computer and on the website.
17. In your terminal, making sure you are in the notebook directory:  
`git add _config.yml` This adds the file you just edited to Git
18. `git commit -m "adding new config file"` Commits your changes and gives a little description of what you changed. You always have to include a commit message, and it is helpful if later you need a specific version of a file and commit messages were well labeled
19. `git push ` This pushes your changes to GitHub.
  _You may have to put in your login information if this is your first push. You can also set this up so you won't have to login again: `git config --global user.name "John Doe"` and `git config --global user.email johndoe@example.com`_
20. There are a few more things to change. On your computer open the _layouts folder and open the default.html file. Edit line 31 to change the link from the Putnam lab link to the Unckless Lab link (http://www.uncklesslab.com/). Edit line 45 to be the link to your GitHub notebook repository (ex. https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/)
21. Save the file and push it to GitHub  
`git add _layouts` Here is a little bit of a shortcut. You can add the whole folder even if you just edited one file. This makes things a lot easier if you edit a lot of posts at once and can commit them all in one line of code. You also don't have to co in and out of directories to add things to git.  
`git commit -m "adding new admin link"`  
`git push origin master`
22. You also have a file called README.md. If you saw at the bottom of the <>Code page on GitHub there is a little window that says Open For Science. This is where Maggie's previous PI has a message about open science and the reasoning behind maintaining an online lab notebook. This is a file that automatically "posts" to the bottom of the repository and usually has a description of it. You should edit it to have your name, a description, and whatever else you want in it. To edit this file open it in your text editor and edit it using Markdown syntax. If you don't know what that is look at these helpful guides: [Markdown file guide](https://guides.github.com/features/mastering-markdown/), [Markdown Cheatsheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet), and [tutorial](https://www.markdowntutorial.com/).
23. Save and add to GitHub in the same way  
`git add README.md`  
`git commit -m "adding new README file"`  
`git push origin master`
24. The last thing you'll probably want to edit before moving on to posts is the about page. It is currently set up as a CV but you can change it to be whatever you'd like. Again this is in Markdown syntax, so take a good look at the example of Maggie's (how the table is written, etc) to see the style.  
`git add about.md`  
`git commit -m "adding new about file"`  
`git push origin master`
25. Now you're almost ready to make a post. However, there are already posts on your page from where you forked your repository. In GitHub you can click on the _posts folder, click on each post, click the small trashcan icon in the upper righthand corner and delete the post. It'll ask you to commit and it looks pretty serious, but as long as you are deleting from your repository (ex. https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/) you are fine. Make sure to not delete the post_template.sh file. As a general rule if you delete all files in a folder it will remove the folder. If you want to do this to quickly remove everything, just remember you'll have to make the _posts folder again. It has to have that exact name so that the files inside will turn into "blog" style posts.
26. Before we can change around things on the computer and push them to GitHub, you have to pull the changes you made in GitHub to your computer. It can be really good practice to always pull from your repository everytime before you want to add something new from your computer. If you are are trying to "merge" two different versions of the repository it can have issues. In the terminal:  
`git pull origin master`
27. To make a new post open a blank file in your text editor and save it to the _posts directory for your notebook. **Important: you have to use the specific Jeykll markdown file naming convention for the posts/blog to build on the site properly: YYYY-MM-DD-Post-Name.md. Also note that you cannot have two different file names be exactly the same except for the date. For example 2022-05-05-DNA-extraction.md and 2022-05-25-DNA-extraction.md will not work. Only one of these notebook pages will post to the website. You have to have unique file names after the date every time.**
28. The first section of your file needs something called a YAML header  
`---`  
`layout: post`  
`title: Your title`  
`category: [ the overarching category of your post ]`  
`tags: [ specfic, tags, if, you, want, them]`  
`---`  
This is what sets up the layout for each post. This jekyll now style blog is set up to track posts with categories and tags. Not all other potential lab notebook sites will have these. If you don't want to use them, just don't add those lines in the YAML header. **Note: do not put a colon (:) in the title, it will make the notebook page not post to the website.**
29. Write your post in Markdown format. Look at [Maggie's Notebook](https://meschedl.github.io/Unckless-Lab-Notebook-Maggie/), and her [notebook pages](https://github.com/meschedl/Unckless-Lab-Notebook-Maggie/tree/master/_posts) in "raw" format by clicking on an md file, and clicking the "raw" button on the top right of the file. This will show you the markdown of the file un-rendered. Looking at the Raw version of any of these files to get the syntax for each element, ex links, images, and tables. Save images into your image folder, you can create other folders for things like pdfs or R scripts. Think of your notebook repository just as your physical one, where you put all your information. Csv files can even be easily added to the repo and looked at online: [example](https://github.com/meschedl/Unckless_Lab_Resources/blob/main/BAC-DiNV/spacer_gen_output.csv). You can look at any of Maggie's notebook repositories [here](https://github.com/meschedl/Unckless-Lab-Notebook-Maggie), [here](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook), and [here](https://github.com/meschedl/MES_Puritz_Lab_Notebook) for ideas. Additionally, here are some suggested guidleines you could follow for writting up a notebook post that is [sample processing](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/Sample-Processing/), or a [reference protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/Protocol-Post-Template/)