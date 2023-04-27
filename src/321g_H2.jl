# Read the contents of the file into a string
json_string = read("/home/arnabbachhar/Intgl/basis/3-21guc_H.json", String)

# Parse the JSON string into a Julia data structure
bs = JSON3.read(json_string)
nested_object_property = bs.elements["1"].electron_shells


