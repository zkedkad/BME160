#!/usr/bin/env python3
# Name: Ziad Kedkad
# Group Members: None
"""
myBio uses classes to show my shortened biography in a line by line output.
The information must be held in constants and not taken in via input.
"""

class Person:
    def __init__(self, name):
        self.myName = name

    def introduce(self):
        print("Hi there, I am ", self.myName)
        name = "Ziad Kedkad"
        username = "dbernick"
        studentType = "undergraduate"


p = Person ('Ziad')
p.introduce()