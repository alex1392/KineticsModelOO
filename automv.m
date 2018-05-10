function automv
  movefile *.mat    mat
  movefile *IF*     IF
  movefile *charts* charts
  movefile *SC*     SC
  %movefile *MS*     MS
  cd charts;  im2mv; movefile *.mp4 ..; cd ..;
  cd IF;      im2mv; movefile *.mp4 ..; cd ..;
  cd SC;      im2mv; movefile *.mp4 ..; cd ..;
  %cd MS;      im2mv; movefile *.mp4 ..; cd ..;
