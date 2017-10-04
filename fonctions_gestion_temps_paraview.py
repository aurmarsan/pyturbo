"""Fonction qui peuvent etre appelees a l'interieur 
d'un programmable filter ou programmable source

remplacement de UPDATE_TIME_STEPS par UPDATE_TIME_STEP
dans la version 3.4
"""
try: from paraview import vtk 
except: import vtk
#__________________________________________________________________________________________
def SetOutputTimesteps(algorithm, timesteps):
    """renseigne la cle qui indique quel instant
    est en sortie du filtre programmable
    
    """
    executive = algorithm.GetExecutive()
    outInfo = executive.GetOutputInformation(0)
    outInfo.Remove(executive.TIME_STEPS())
    for timestep in timesteps:
        outInfo.Append(executive.TIME_STEPS(), timestep)
        outInfo.Remove(executive.TIME_RANGE())
        outInfo.Append(executive.TIME_RANGE(), timesteps[0])
        outInfo.Append(executive.TIME_RANGE(), timesteps[-1])
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def GetUpdateTimestep(algorithm):
    """retourne l'instant demande par le 
    pipeline paraview
    
    """
    executive = algorithm.GetExecutive()
    outInfo = executive.GetOutputInformation(0)
    if not outInfo.Has(executive.UPDATE_TIME_STEP()):
        return None
    return outInfo.Get(executive.UPDATE_TIME_STEP())
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def get_temps_animation(filtre_programmable):
    """fonction qui appelle GetUpdateTimestep
    
    retourne la valeur de temps demande par le pipeline Paraview
    a utiliser dans la partie SCRIPT du filtre programmable
    
    """
    temps_demande = GetUpdateTimestep(filtre_programmable)
    return temps_demande
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def set_sortie_animable(filtre_programmable, liste_temps=None):
    """fonction qui indique que le filtre programmable produit une sortie
    dependante du temps demande par le pipeline paraview
    
    doit etre mis dans la partie REQUESTINFORMATION du filtre programmable
    sinon, pas de mise a jour automatique de la source
    
    There's one caveat. Currently ParaView client does not notice the timesteps produced 
    by the source and hence the GUI will not update the animation time ranges automatically, .
    as is the case for standard reader. 
    http://paraview.org/Wiki/Here_are_some_more_examples_of_simple_ParaView_3_python_filters.#Producing_Data_with_Timesteps_.28Source.29
    
    on met donc (0, 10, 20, 30) ... 
    
    liste_temps peut toutefois etre indique dans le cas ou l'on souhaite ensuite suavegarder 
    les donnees sous forme de vtm, c'ets cette liste qui est utilisee"
    """
    if liste_temps is not None:
        SetOutputTimesteps(filtre_programmable, tuple(liste_temps))
    else:
        SetOutputTimesteps(filtre_programmable, (0, 10, 20, 30))
    return 0
#__________________________________________________________________________________________

#__________________________________________________________________________________________
def set_output(self, vtkDataObject):
    """fonction pour indiquer la sortie d'un filtre ou source programmable
    """
    #la sortie polydata, entre autres, n'est pas bien gere par le filtre programmable
    #on transforme alors dans ce cas en multibloc d'un seul bloc, qui marche bien. 
    if not isinstance(vtkDataObject, vtk.vtkMultiBlockDataSet):
        m = vtk.vtkMultiBlockDataSet()
        m.SetBlock(0, vtkDataObject)
        vtkDataObject = m
    #on met en sortie
    self.GetOutput().ShallowCopy(vtkDataObject)
#__________________________________________________________________________________________
