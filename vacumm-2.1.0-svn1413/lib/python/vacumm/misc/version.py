"""Get version of misc tools using svn revision number of directory"""
import subprocess, os
target_path = os.path.dirname(os.path.abspath(__file__))
version = subprocess.Popen(["svnversion", target_path], stdout=subprocess.PIPE).communicate()[0].split(':')[-1].split('\n')[0]

