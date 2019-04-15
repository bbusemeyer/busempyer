import numpy as np

class Scheduler:
  def __init__(self,first_date,presenters,avoid=[]):
    # Options.
    self.avoid=avoid
    self.next_date=first_date
    self.presenters=presenters

    # Hard code.
    self.months = ["Jan.","Feb.","March","April","May","June","July","August","Sept.","Oct.","Nov.","Dec."]
    self.day_month = [31, 28,    31,     30,     31,   30,    31,    31,      30,     31,    30,    31]
    self.midx=self.months.index(first_date[0])
    self.dates=[]

  def _step_date(self):
    self.next_date[1] += 7
    if self.next_date[1] > self.day_month[self.midx]:
      self.next_date[1] -= self.day_month[self.midx]
      self.midx += 1
      self.next_date[0] = self.months[self.midx]


  def find_dates(self):
    for person in self.presenters:
      self.dates.append("{0:<6} {1:>2}".format(*self.next_date))
      self._step_date()
      while self.next_date in self.avoid:
        self._step_date()

  def shuffle(self):
    """ Exports a suggested schedule. New for every call."""
    np.random.shuffle(self.presenters)

  def export(self):
    for date,presenter in zip(self.dates,self.presenters):
      print("{}   {}".format(date,presenter))

if __name__=='__main__':
  avoid = [["Jan.", 21]]
  presenters = [
      'Will','Brian','Shivesh','João','Kiel','XioLei','Yueqing','Mick','Chong','Alex','Paul'
    ]
  first_date = ["Jan.",14]
  sched=Scheduler(first_date,presenters,avoid)
  sched.find_dates()
  sched.export()
