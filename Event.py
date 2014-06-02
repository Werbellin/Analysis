class Event :
    def __init__(self, leaf) :
        self.data = leaf
        self.Cuts = {'start':True}
        self.CutHistory = ""
