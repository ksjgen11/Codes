class Point:
    x = 0
    y = 0
    
    def __init__(self, x=0, y=0):       # Constructor with default init values
        self.x = x
        self.y = y
    
    def translate(self, dx, dy):
        self.x += dx
        self.y += dy
        
    def __str__(self):                  # Now objects can respond to "string requests"
        return 'Point (%.1f, %.1f)' % (self.x, self.y)
    
    def __add__(self, other):           # Over-riding +
        p = Point(self.x, self.y)
        p.translate(other.x, other.y)
        return p
    
origo = Point()


    