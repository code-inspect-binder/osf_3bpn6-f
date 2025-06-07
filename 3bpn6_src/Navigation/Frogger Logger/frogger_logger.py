# coding: UTF-8

'''
Frog's POsition & Direction Coder v.0.1
jinook.oh@univie.ac.at
Cognitive Biology Dept., University of Vienna
- 2017.1

----------------------------------------------------------------------
Copyright (C) 2016 Jinook Oh, W. Tecumseh Fitch 
- Contact: jinook.oh@univie.ac.at, tecumseh.fitch@univie.ac.at

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
----------------------------------------------------------------------

Usage
1) Determine filename and write it down on the text-control on the 
upper left corner.
2) Click(press) 'Start' button
3) Click(press) to determine frog's position and then click again
to determine frog's direction
4) Can go back to a previous position to delete or adjust by using
buttons at the bottom
5) Click 'Stop' button to save and stop coding.
6) Click 'Open' button to open previous result CSV file and adjust
it. (Timestamp will not be changed)
7) Click 'Quit' to quit the app
'''

from os import getcwd, path
from sys import argv
from copy import copy
from time import time, sleep
from datetime import timedelta
from glob import glob

import wx

from misc_funcs import GNU_notice, get_time_stamp, writeFile, show_msg, load_img, calc_line_angle, calc_pt_w_angle_n_dist

# ======================================================

class FPODICFrame(wx.Frame):
    def __init__(self):
        self.w_size = (500, 700)
        wx.Frame.__init__(self, None, -1, 'Frog POsition & DIrection Coder', size=self.w_size, style=wx.DEFAULT_FRAME_STYLE) # init frame
        self.SetPosition( (0, 20) )
        self.Show(True)
        self.panel = wx.Panel(self, pos=(0,0), size=self.w_size)
        self.panel.SetBackgroundColour('#000000')

        ### init variables
        self.program_start_time = time()
        self.session_start_time = -1
        self.oData = [] # output data
        self.di = 0 # current data index
        self.fPath = '' # result file path
        self.timer_run = None # timer for running analysis
        self.cStat = '' # to indicate a mouse click is for positioning or direction.
        self.flag_new = True # indicating whether it's a new session or opened session (opened session timestampe will not be recorded)
        self.flag_adj_grid = False # grid adjusting
        self.m_press_pt = (-1,-1)
        self.num_of_circles = 4
        self.angle_of_guide_lines = 20
        
        ### upper user interface setup
        posX = 5
        posY = 5
        self.txt_name = wx.TextCtrl(self.panel, -1, value='Enter filename here', pos=(posX,posY), size=(150,-1)) # output file name
        self.txt_name.Bind(wx.EVT_LEFT_UP, self.onTextCtrl)
        posX += self.txt_name.GetSize()[0]
        posY = 2
        self.btn_start = wx.Button(self.panel, -1, name='start', label='Start', pos=(posX,posY), size=(50, -1))
        self.btn_start.Bind(wx.EVT_LEFT_UP, self.onStartStopAnalyze)
        posX += self.btn_start.GetSize()[0] + 30
        self.btn_grid = wx.Button(self.panel, -1, name='grid', label='Adj. Grid OFF', pos=(posX,posY), size=(100, -1))
        self.btn_grid.Bind(wx.EVT_LEFT_UP, self.onAdjustGrid)
        posX = self.w_size[0]-120 
        posY = 2
        self.btn_open = wx.Button(self.panel, -1, name='open', label='Open', pos=(posX,posY), size=(50, -1))
        self.btn_open.Bind(wx.EVT_LEFT_UP, self.onStartStopAnalyze)
        posX += self.btn_open.GetSize()[0]+5
        self.btn_quit = wx.Button(self.panel, -1, name='quit', label='Quit', pos=(posX,posY), size=(50, -1))
        self.btn_quit.Bind(wx.EVT_LEFT_UP, self.onClose)
        posX = 5
        posY += self.txt_name.GetSize()[1]+10
        self.font = wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL)
        self.sTxt_pr_time = wx.StaticText(self.panel, -1, label='0:00:00', pos=(posX, posY)) # elapsed time since program starts
        self.sTxt_pr_time.SetFont(self.font)
        posX = self.sTxt_pr_time.GetPosition()[0] + self.sTxt_pr_time.GetSize()[0] + 10
        _stxt = wx.StaticText(self.panel, -1, label='(program)', pos=(posX, posY))
        _stxt.SetForegroundColour('#CCCCCC')
        _stxt.SetFont(self.font)
        self.sTxt_pr_time.SetBackgroundColour('#000000')
        self.sTxt_pr_time.SetForegroundColour('#00FF00')
        posX += _stxt.GetSize()[0] + 30
        self.sTxt_s_time = wx.StaticText(self.panel, -1, label='0:00:00', pos=(posX, posY)) # elapsed time since session starts
        self.sTxt_s_time.SetFont(self.font)
        posX = self.sTxt_s_time.GetPosition()[0] + self.sTxt_s_time.GetSize()[0] + 10
        _stxt = wx.StaticText(self.panel, -1, label='(session)', pos=(posX, posY))
        _stxt.SetForegroundColour('#CCCCCC')
        _stxt.SetFont(self.font)
        self.sTxt_s_time.SetBackgroundColour('#000000')
        self.sTxt_s_time.SetForegroundColour('#CCCCFF')

        ### data coding panel 
        posY += self.sTxt_s_time.GetSize()[1]+10
        self.cp_pos = (5, posY)
        self.cp_sz = (self.w_size[0]-15, self.w_size[0]-15)
        self.cp = wx.Panel(self.panel, pos=self.cp_pos, size=self.cp_sz)
        self.cp.Bind(wx.EVT_PAINT, self.onPaint)
        self.cp.Bind(wx.EVT_LEFT_DOWN, self.onMouseLeftDown)
        self.cp.Bind(wx.EVT_MOTION, self.onMouseMove)
        self.cp.Bind(wx.EVT_LEFT_UP, self.onMouseLeftUp)

        ### navigation and delete button
        _w = 50
        posX = (self.w_size[0])/2 - int(_w*2.5+30)
        posY += self.cp_sz[1]+5
        self.sTxt_info = wx.StaticText(self.panel, -1, label='[Idx. 00] / timestamp 00:00:00  /  x 0.000 / y 0.000  /  dir 000', pos=(posX, posY)) # several information display
        self.sTxt_info.SetPosition((self.w_size[0]/2-self.sTxt_info.GetSize()[0]/2, posY))
        self.sTxt_info.SetForegroundColour('#CCCCCC')
        posY += self.sTxt_info.GetSize()[1]+10 
        self.btn_prev2 = wx.Button(self.panel, -1, name='prev2', label='<<', pos=(posX,posY), size=(_w, -1))
        self.btn_prev2.Bind(wx.EVT_LEFT_UP, self.onNavigate)
        posX += self.btn_prev2.GetSize()[0]+10
        self.btn_prev = wx.Button(self.panel, -1, name='prev', label='<', pos=(posX,posY), size=(_w, -1))
        self.btn_prev.Bind(wx.EVT_LEFT_UP, self.onNavigate)
        posX += self.btn_prev.GetSize()[0]+10
        self.btn_next = wx.Button(self.panel, -1, name='next', label='>', pos=(posX,posY), size=(_w, -1))
        self.btn_next.Bind(wx.EVT_LEFT_UP, self.onNavigate)
        posX += self.btn_next.GetSize()[0]+10
        self.btn_next2 = wx.Button(self.panel, -1, name='next2', label='>>', pos=(posX,posY), size=(_w, -1))
        self.btn_next2.Bind(wx.EVT_LEFT_UP, self.onNavigate)
        posX += self.btn_next2.GetSize()[0]+20
        self.btn_del = wx.Button(self.panel, -1, label='Del.', pos=(posX,posY), size=(_w, -1))
        self.btn_del.Bind(wx.EVT_LEFT_UP, self.onRemove)
       
        statbar = wx.StatusBar(self, -1)
        self.SetStatusBar(statbar)

        ### set timer for updating the current running time
        self.timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.onTimer, self.timer)
        self.timer.Start(1000)

        self.Bind( wx.EVT_CLOSE, self.onClose )

    # --------------------------------------------------       
    
    def onPaint(self, event):
        dc = wx.PaintDC(event.GetEventObject())
        self.draw_graph(dc) 
        event.Skip()
    
    # --------------------------------------------------       
    
    def draw_graph(self, dc):
        dc.SetBackground(wx.Brush(wx.Colour(50,50,50)))
        dc.Clear()
        
        dc.SetPen(wx.Pen(wx.Colour(100,100,100), 1))
        ### draw direction guide lines
        for angle in range(0,360,self.angle_of_guide_lines):
            if angle % 90 == 0: continue
            pt = list( calc_pt_w_angle_n_dist(angle, self.cp_sz[0]/2) )
            pt[0] += self.cp_sz[0]/2
            pt[1] = self.cp_sz[1]/2 - pt[1]
            dc.DrawLine(self.cp_sz[0]/2, self.cp_sz[1]/2, pt[0], pt[1])
        dc.SetBrush(wx.Brush(wx.Colour(0,0,0), wx.TRANSPARENT))
        ### draw guide circles (for distance)
        _r = self.cp_sz[0]/2/self.num_of_circles
        for r in range(_r, self.cp_sz[0]/2+1, _r):
            dc.DrawCircle(self.cp_sz[0]/2, self.cp_sz[1]/2, r)
        ### draw major direction guide lines
        dc.SetPen(wx.Pen(wx.Colour(200,200,200), 1))
        dc.DrawLine(self.cp_sz[0]/2, 0, self.cp_sz[0]/2, self.cp_sz[1]) 
        dc.DrawLine(0, self.cp_sz[1]/2, self.cp_sz[0], self.cp_sz[1]/2) 
        dc.DrawLine(0, 0, self.cp_sz[0], self.cp_sz[1]) 
        dc.DrawLine(self.cp_sz[0], 0, 0, self.cp_sz[1])

        if len(self.oData) > 0:
            ### draw past data points 
            c = 30
            inc = (255-30)/len(self.oData)
            for i in range(len(self.oData)):
                if i==len(self.oData)-1 and len(self.oData[i])==1: break # the last one without direction. it'll be drawn in below 
                dc.SetPen(wx.Pen(wx.Colour(200,200,200), 1))
                dc.SetBrush(wx.Brush(wx.Colour(c,c,c)))
                _d = self.oData[i]
                c += inc 
                if i == self.di: # current one
                    ### draw extra circle to indicate the current data 
                    dc.SetPen(wx.Pen(wx.Colour(255,255,0), 1))
                    dc.SetBrush(wx.Brush(wx.Colour(c,c,c), wx.TRANSPARENT))
                    dc.DrawCircle(_d['pos'][0], _d['pos'][1], 7)
                    dc.SetPen(wx.Pen(wx.Colour(200,200,200), 1))
                    dc.SetBrush(wx.Brush(wx.Colour(c,c,c)))
                dc.DrawCircle(_d['pos'][0], _d['pos'][1], 5)
                pt = list( calc_pt_w_angle_n_dist(_d['dir'], 10) )
                pt[0] += _d['pos'][0]
                pt[1] = _d['pos'][1]-pt[1]
                dc.SetPen(wx.Pen(wx.Colour(50,255,50), 2))
                dc.DrawLine(_d['pos'][0], _d['pos'][1], pt[0], pt[1])
                if i > 0:
                    dc.SetPen(wx.Pen(wx.Colour(c,c,0), 1))
                    _pd = self.oData[i-1]
                    dc.DrawLine(_pd['pos'][0], _pd['pos'][1], _d['pos'][0], _d['pos'][1])
                
            ### draw guide lines to determine direction 
            if self.cStat == 'to_click_dir':
                cPt = self.oData[self.di]['pos']
                dc.SetPen(wx.Pen(wx.Colour(0,200,0), 1))
                for angle in range(0,360,20):
                    pt = list( calc_pt_w_angle_n_dist(angle, 100) )
                    pt[0] += cPt[0] 
                    pt[1] = cPt[1] - pt[1]
                    dc.DrawLine(cPt[0], cPt[1], pt[0], pt[1])
                dc.SetPen(wx.Pen(wx.Colour(100,100,100), 1))
                dc.SetBrush(wx.Brush(wx.Colour(255,150,0)))
                dc.DrawCircle(cPt[0], cPt[1], 10)
       
            if self.oData[self.di].has_key('timestamp') == True:
                ### display some info of the current data
                _pos = self.calc_nav_pos(self.di)
                self.sTxt_info.SetLabel('[Idx. %i] / timestamp %s  /  x %.3f / y %.3f  /  dir %i'%(self.di, 
                                                                                                 self.oData[self.di]['timestamp'], 
                                                                                                 _pos[0], 
                                                                                                 _pos[1], 
                                                                                                 self.oData[self.di]['dir']))
    
    # --------------------------------------------------       
    
    def calc_nav_pos(self, di):
        ''' converts a data point's position (pixel position on screen)
        to a navigation position coordinates; (0,0) is the center & 
        right/up most pixel is 1.0 & left/down most pixel is -1.0
        '''
        _x = (self.oData[di]['pos'][0]-self.cp_sz[0]/2) / float(self.cp_sz[0]/2)
        _y = (self.cp_sz[1]/2-self.oData[di]['pos'][1]) / float(self.cp_sz[1]/2)
        return (_x, _y)

    # --------------------------------------------------       
    
    def onTextCtrl(self, event):
        self.txt_name.SetSelection(-1, -1)
        event.Skip()
    
    # --------------------------------------------------       
    
    def onAdjustGrid(self, event):
        self.m_press_pt = (-1,-1)
        if self.flag_adj_grid == False:
            self.flag_adj_grid = True
            self.btn_grid.SetLabel('Adj.Grid ON')
        elif self.flag_adj_grid == True:
            self.flag_adj_grid = False 
            self.btn_grid.SetLabel('Adj.Grid OFF')
        event.Skip()
    
    # --------------------------------------------------       
    
    def onMouseLeftDown(self, event):
        if self.flag_adj_grid == True:
            self.m_press_pt = event.GetPosition()
        event.Skip()
    
    # --------------------------------------------------       
    
    def onMouseMove(self, event):
        ### currently adjusting number of grid lines.
        if self.flag_adj_grid==False or self.m_press_pt==(-1,-1): return 
        mp = event.GetPosition()
        r_len = self.cp_sz[0]/2 # total ref. length
        _w = mp[0]-self.m_press_pt[0] # horizontal move from mouse button press point
        _h = self.m_press_pt[1]-mp[1] # vertical move
        u_len = r_len / 10 # unit length to determine the number of circles, which could be 10 at most
        self.num_of_circles = min(max(2, _w/u_len), 10)
        u_len = r_len / 2 # unit length to determine the angle between two direction guide lines, which could be 10 or 20. 
        self.angle_of_guide_lines = min(max(10, (2-_h/u_len)*10), 20)
        self.cp.Refresh()
        event.Skip()
    
    # --------------------------------------------------       
    
    def onMouseLeftUp(self, event):
        if self.flag_adj_grid == True:
            self.m_press_pt = (-1,-1)
            return
        if self.fPath == '': return
        mp = event.GetPosition()
        if self.cStat == 'to_click_pos':
            ### determine position
            if self.di == len(self.oData)-1:
                self.oData.append({})
                self.di += 1
            self.oData[self.di]['pos'] = mp
            self.cStat = 'to_click_dir'
        elif self.cStat == 'to_click_dir': # clicked to determine direction
            ### determine direction
            self.oData[self.di]['dir'] = calc_line_angle(self.oData[self.di]['pos'], mp)
            self.cStat = 'to_click_pos'
            if (self.flag_new==True and self.di==len(self.oData)-1) or (self.flag_new==False and self.di==len(self.oData)-1 and self.oData[self.di].has_key('timestamp')==False): self.oData[self.di]['timestamp'] = self.sTxt_s_time.GetLabel()
        self.cp.Refresh()
        event.Skip()
    
    #------------------------------------------------
 
    def onNavigate(self, event):
        ''' select a data point
        '''
        if self.fPath =='' or len(self.oData)==0: return
        en = event.GetEventObject().GetName()
        if en == 'prev2': self.di = max(0, self.di-10)
        elif en == 'prev': self.di = max(0, self.di-1)
        elif en == 'next': self.di = min(self.di+1, len(self.oData)-1)
        elif en == 'next2': self.di = min(self.di+10, len(self.oData)-1)
        self.cp.Refresh()
        event.Skip()

    #------------------------------------------------
    
    def onRemove(self, event):
        ''' remove a data point
        '''
        self.oData.pop(self.di)
        if self.di >= len(self.oData): self.di=len(self.oData)-1
        self.cp.Refresh()
        event.Skip()
    
    #------------------------------------------------
    
    def onStartStopAnalyze(self, event):
        en = event.GetEventObject().GetName()
        if self.session_start_time == -1: # not in analysis session. start a session
            if en == 'start':
                fn = self.txt_name.GetValue().strip()
                if fn=='Enter filename here' or fn=='':
                    show_msg('Please enter proper filename to save results')
                    return
                self.fPath = path.join(CWD, fn+'.csv')
                self.di = -1
                self.oData = []
                self.flag_new = True
            elif en == 'open':
                dlg = wx.FileDialog(self, "Select CSV result file to open.", CWD, wildcard='(*.csv)|*.csv', style=wx.FD_DEFAULT_STYLE|wx.FD_FILE_MUST_EXIST)
                if dlg.ShowModal() == wx.ID_CANCEL: return
                self.fPath = dlg.GetPath()    
                self.di = 0
                self.oData = []
                ### read data from the selected file
                fh = open(self.fPath, 'r')
                lines = fh.readlines()
                for line in lines:
                    items = [x.strip() for x in line.split(',')]
                    if len(items) < 5: continue
                    idx = None
                    try: idx = int(items[0])
                    except: pass
                    if idx != None:
                        self.oData.append( dict(timestamp=items[1], 
                                                pos=(int(items[2]),int(items[3])), 
                                                dir=int(items[4])) )
                fh.close()
                self.txt_name.SetValue( path.basename(self.fPath).replace('.csv','') )
                self.flag_new = False
            self.session_start_time = time()
            self.btn_start.SetLabel('Stop')
            self.cStat = 'to_click_pos'
        else: # in session. stop it.
            result = show_msg(msg='Save current data?', cancel_btn = True)
            if result == True: self.onSave()
            self.session_start_time = -1
            self.sTxt_s_time.SetLabel('0:00:00')
            self.btn_start.SetLabel('Start')
            self.txt_name.SetValue('Enter filename here')
            self.di = -1
            self.oData = []
        self.cp.Refresh()
        event.Skip()
    
    # --------------------------------------------------       
    
    def onSave(self):
        ### save CSV result file
        fn = self.txt_name.GetLabel() + '.csv'
        fp_n = path.basename(self.fPath)
        if fp_n != fn: self.fPath.replace(fp_n, fn)
        fh = open(self.fPath, 'w')
        fh.write('data-index, Timestamp, posX, posY, Direction\n')
        for di in range(len(self.oData)):
            d_ = self.oData[di]
            line = '%i, %s, %s, %s, %s\n'%(di,
                                           d_['timestamp'],
                                           str(d_['pos'][0]), 
                                           str(d_['pos'][1]), 
                                           str(d_['dir'])
                                           )
            fh.write(line)
        #fh.write('------------------------------------------------------------------\n')
        fh.close()

        ### save image
        bmp = wx.EmptyBitmap(self.cp_sz[0], self.cp_sz[1], depth=-1)
        memdc = wx.MemoryDC()
        memdc.SelectObject(bmp)
        self.draw_graph(memdc)
        memdc.SelectObject(wx.NullBitmap)
        img = bmp.ConvertToImage()
        fp = self.fPath.replace('.csv', '.png')
        img.SaveFile(fp, wx.BITMAP_TYPE_PNG)

        msg = 'Saved.\n'
        chr_num = 50 # characters in one line
        if len(self.fPath) > chr_num:
            for i in range(len(self.fPath)/chr_num):
                msg += '%s\n'%(self.fPath[chr_num*i:chr_num*(i+1)])
            msg += '%s\n'%(self.fPath[chr_num*(i+1):])
        else:
            msg += '%s'%(self.fPath)
        show_msg(msg, size=(400,200))

    # --------------------------------------------------       

    def onTimer(self, event):
        ''' Main timer 
        updating running time on the main window
        '''
        ### update several running time
        e_time = time() - self.program_start_time
        self.sTxt_pr_time.SetLabel( str(timedelta(seconds=e_time)).split('.')[0] )
        if self.session_start_time != -1:
            e_time = time() - self.session_start_time
            self.sTxt_s_time.SetLabel( str(timedelta(seconds=e_time)).split('.')[0] )

    # --------------------------------------------------

    def show_msg_in_statbar(self, msg, time=5000):
        self.SetStatusText(msg)
        wx.FutureCall(time, self.SetStatusText, "") # delete it after a while

    # --------------------------------------------------

    def onClose(self, event):
        self.timer.Stop()
        result = True
        if self.session_start_time != -1: # session is running
            result = show_msg(msg='Session is not stopped..\nUnsaved data will be lost. (Stop analysis to save.)\nOkay to exit anyway?', cancel_btn = True)
        if result == True: wx.FutureCall(100, self.Destroy)

# ======================================================

class FPODICApp(wx.App):
    def OnInit(self):
        self.frame = FPODICFrame()
        self.frame.Show()
        self.SetTopWindow(self.frame)
        return True

# ======================================================

if __name__ == '__main__':
    if len(argv) > 1:
        if argv[1] == '-w': GNU_notice(1)
        elif argv[1] == '-c': GNU_notice(2)
    else:
        GNU_notice(0)
        CWD = getcwd()
        app = FPODICApp(redirect = False)
        app.MainLoop()




