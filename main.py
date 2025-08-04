# Script to automatically run the three programs

import subprocess

subprocess.run(["python", "radiator.py"])
subprocess.run(["python", "hybrid.py"])
subprocess.run(["python", "pcm.py"])
subprocess.run(["python", "pcmProperty.py"])
subprocess.run(["python", "post.py"])