class A(object):
    def __init__(self):
        self.x = 'Hello'

    def method_a(self, foo):
        print(self.x + ' ' + foo)

me = A()
print(me.x)
me.method_a("Fred")

