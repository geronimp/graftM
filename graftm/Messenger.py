#!/usr/bin/env python
from datetime import datetime

class Messenger:
    
    def __init__(self):
        self.time = datetime.now().strftime('%H:%M:%S')
        self.date = datetime.now().strftime('%d-%m-%Y')
        
    def message(self, message):
        print '\t%s  [%s]: %s' % (self.date, self.time, message)
    
    def error_message(self, message):
        print '\t%s  [%s]: == ERROR == %s' % (self.date, self.time, message)
        
    def header(self, message):
        print '\n[%s]: %s' % (self.time, message)