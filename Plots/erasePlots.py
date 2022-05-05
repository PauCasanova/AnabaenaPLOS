import os

folder=os.getcwd()
test = os.listdir(folder)
for images in test:
	if images.endswith(".png"):
		os.remove(os.path.join(folder, images))
