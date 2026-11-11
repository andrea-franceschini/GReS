function chronos(val)

switch lower(val)
  case 'on'
    
  case 'off'
  case 'force'
  otherwise
    error('Unknown option for chronos setup')
end

end

