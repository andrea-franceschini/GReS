# Controlla esistenza Output/input_ns3d e le crea nel caso non ci siano
[ -d Output/input_ns3d ] || mkdir -p Output/input_ns3d

# Avvio simulazione
nohup ./biogeo > screen &
