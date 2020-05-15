from __future__ import print_function
import cobra
import cobra.test


# "ecoli" and "salmonella" are also valid arguments
model = cobra.test.create_test_model("textbook")

print(len(model.reactions))
print(len(model.metabolites))
print(len(model.genes))

# Get out one reaction
model.reactions[29]

# Use dictionary to get item
model.metabolites.get_by_id("atp_c")