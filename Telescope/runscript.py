# Script to build and execute c++ tracer code, and then postprocess in python.
# Usage:
#choose what you want to build and run, plus add your program option

# python runscript.py <makerule> <executable> <program_option> <postprocess_option>'

import sys, os, getopt
import subprocess
#from subprocess import run,PIPE

def main(argv): #takes list of arguments

#If all python argv are passed accordingly.
    if len(argv)==4: 
        make_rule    = argv[0]
        exe_file     = argv[1]
        exe_argv     = argv[2]
        post_process = argv[3]
        try:
            print("\nBuilding Program..\n")
            subprocess.run(['make', make_rule], check=True)
        except subprocess.CalledProcessError:
            print("\nError in building\n")
            help()
        try:
            print("\nRunning Program...\n")
            subprocess.run(["./"+exe_file+" "+exe_argv], shell=True, check = True)
        except subprocess.CalledProcessError:
            print("\nError in running c++ executable\n")
            help()

        if(post_process == 1):
            try:
                print("\nPost Processing and Plots...\n")    
                subprocess.run(["python3","Post_processing_and_Plots.py"], check = True)
            except subprocess.CalledProcessError:
                print("\nError in Post Processing\n")
                help()
        else:
            print("\nNo further processing, no new plots\n")
#Clean all
    elif (len(argv)==1 and argv[0]=="reset"):
        os.system("make allclean")

#If no valid argcount.
    else: 
        help()  

        
    
    


def help():
    print("Try running as:")
    print("python <script> <make_rule> <exe_file> <exe_argv> <post_process>")
    print("make_rule    -> may produce exe with same name.")
    print("exe_file     -> executable")
    print("exe_argv     -> argument variable for simulation type.")
    print("post_process -> post processing option(1=yes, 0=no)")
    print("To reset run as:")
    print("python 'reset' ")



#script
if __name__=='__main__':
    main(sys.argv[1:])









    """
    def main(): #takes list of arguments

        target_dict = {0:"default", 1:"tracer", 2:"ray_int", 3:"allclean", 4:"clean"}
        target = input("Choose building process:\n0.default\n1.tracer\n2.ray_int\n3.allclean\n4.clean\n99.help")
        if( target in target in target_dict):
            os.system("make " + target)
        else:
            help_make()
        

        print("Running Program...")
        os.system("./"+exe_file +" "+ exe_argv)
        

        
    
    


def help_make():
    print("0.default: Will build everything\n 1.tracer)



#script
if __name__=='__main__':
    main()
"""