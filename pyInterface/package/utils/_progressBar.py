
import sys

class progressBar:

	minimum = 0.0
	maximum = 100.0
	stars = 0
	fp = None
	full = False
	effRange = 0.

	def __init__(self, minimum = 0, maximum = 100, fp=sys.stdout):
		self.fp = fp
		self.reset(minimum, maximum)

	def reset(self, minimum = 0, maximum = 100):
		self.minimum = float(minimum)
		self.maximum = float(maximum)
		self.effRange = 51. / (self.maximum - self.minimum)
		self.minimum = self.minimum - (0.5 / self.effRange)
		self.stars = 0
		self.full = False

	def start(self):
		self.fp.write("0%   10   20   30   40   50   60   70   80   90   100%\n")
		self.fp.write("|----|----|----|----|----|----|----|----|----|----|\n")

	def cancel(self):
		self.fp.write("<!\n")

	def update(self, prog):
		if not self.full:
			stars = int((prog - self.minimum) * self.effRange)
			diff = stars - self.stars
			if diff > 0:
				self.fp.write(diff * "*")
				self.fp.flush()
				self.stars += diff
				if self.stars == 51:
					self.fp.write('\n')
					self.full = True

