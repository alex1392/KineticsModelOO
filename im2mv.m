function im2mv
files = dir('*.jpg');
p = 1;
while numel(files) > 0
  i = 1;
  frame = struct('cdata',[],'colormap',[]);
  while i < 1000 && i <= numel(files)
    try
      fName = files(i).name;
      s = [strfind(fName,'C'),strfind(fName,'F'),strfind(fName,'s')];
      f = strfind(fName,'.');
      ID = str2num(fName(s(end)+1:f(end)-1));
      frame(ID).cdata = imread(files(i).name);
      i = i + 1;
    catch
      files(i) = [];
    end
  end
  fName = files(i-1).name;
  s = strfind(fName,num2str(ID));
  fName = fName(1:s-1);
  if i <  numel(files)
    fName = [fName,'-part',num2str(p)];    
  end
  movie = VideoWriter(fName,'MPEG-4');
  movie.Quality = 100;
  movie.FrameRate = 10;
  open(movie);
  writeVideo(movie, frame);
  close(movie);
  p = p + 1;
  files(1:i-1) = [];
end

