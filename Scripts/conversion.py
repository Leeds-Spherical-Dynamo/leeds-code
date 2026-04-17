"""Short script for converting between leeds dynamo codes.

This script will take in parameters used for a dynamo simulation
in the old leeds dynamo code, and output the equivalent parameters
for the new dynamo code, and vice-versa.
"""
import signal
import sys

def handler(signum, frame):
    res = input("Ctrl-c was pressed. Do you really want to exit? y/n ")
    if res == 'y':
        exit(1)

#signal.signal(signal.SIGINT, handler)

def getInput():
    end=False
    while not end:
        try:
            conversion_type = input("Please enter 1 to convert from Old Dynamo code to new,\n or 2 to convert from New Dynamo code to old:")
            conversion_type=int(conversion_type)
            if conversion_type==1 or conversion_type==2:
                end=True
            else:
                print("please enter 1 or 2")
        except KeyboardInterrupt:
            print("\nExiting...\n")
            sys.exit(0)
        except:
            print("please try again")
    return conversion_type

def getParam(param_name):
    end=False
    while not end:
        try:
            param=float(input("Enter {}:".format(param_name)))
            if param > 0:
                end=True
            else:
                print("{}<0, please try again".format(param_name))
        except KeyboardInterrupt:
            print("\nExiting...\n")
            sys.exit(0)
        except:
            print("input error, try again")
    return param

conversion_type=getInput()


if conversion_type==1:
    #converting from old to new
    print("---Converting Old Dynamo code parameters to new dynamo code parameters---")
    #raw=float(args['--ra']); ekw=float(args['--ek']); q=float(args['--q']); ro=float(args['--ro']);
    raw=getParam("Ra")
    ekw=getParam("Ek")
    q=getParam("q")
    ro=getParam("Ro")

    raj = raw/ekw
    ekj = 2* ekw
    pm = ekw/ro
    pr = ekw/(q*ro)
    print("\n--------------------------\n")
    print(" Ra:{:.4e}\n Ek:{:.4e}\n Pm:{:.4e}\n Pr:{:.4e}\n".format(raj,ekj,pm,pr))
else:
    #converting from new to old
    print("---Converting new Dynamo code parameters to old dynamo code parameters---")
    #raj=float(args['--ra']); ekj=float(args['--ek']); pr=float(args['--pr']); pm=float(args['--pm']);
    raj=getParam("Ra")
    ekj=getParam("Ek")
    pr=getParam("Pr")
    pm=getParam("Pm")
    raw = (raj*ekj)/2
    ekw = ekj/2
    q = pm/pr
    ro= ekw/(q*pr)
    print("\n--------------------------\n")
    print(" Ra:{:.4e}\n Ek:{:.4e}\n q:{:.4e} \n Ro:{:.4e}\n".format(raw,ekw,q,ro))
