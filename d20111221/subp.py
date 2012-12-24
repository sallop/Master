#!/usr/bin/python
#-*- encoding: utf-8 -*-
import subprocess
import os

# subprocess.Popen(args, bufsize=0, executable=None,
#                  stdin=None, stdout=None, stderr=None,
#                  preexec_fn=None, close_fds=False,
#                  shell=False, cwd=None, env=None,
#                  universal_newlines=False, startupinfo=None, creationflags=0)
# output=`dmesg | grep hda`
# ==>
# p1 = Popen(["dmesg"], stdout=PIPE)
# p2 = Popen(["grep", "hda"], stdin=p1.stdout, stdout=PIPE)
# output = p2.communicate()[0]

# proc = subprocess.Popen(['ruby','read_ublas.rb',
#                          'L3x3-mat.dat',
#                          'L3x3-t_gibbs.dat'] ,
#                         stdin =subprocess.PIPE,
#                         stdout=subprocess.PIPE,)

# while True:
#     next = sys.stdin.readline()
#     if not next: break
#     sys.stdout.write(next)
#     sys.stdout.flush()
# sys.stderr.write('repeater.py: existing\n')
# sys.stderr.flush()


proc = subprocess.Popen(['python','read_ublas.py','dir-dat/L3x3-mat.dat'],
                        stdin =subprocess.PIPE,
                        stdout=subprocess.PIPE,)

(stdout, stdin) = (proc.stdout, proc.stdin)

print "-"*80
# while True:
#     line = stdout.readline()
#     if not line:
#         break
#     print line.strip()

for line in stdout:
    print line.strip()

print "-"*80

ret = proc.wait()
print "Return code:%d"%ret


# output = proc.communicate()[0]
# lines = output.split()

# for line in lines:
#     print line

