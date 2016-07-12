""" 
PROXIMITY RIP
Description: Used for formatting and summarizing Tongjia's Proximity Analysis
data output that has been processed by Derek's Proximity Filter Program. 

NOTE: DOES NOT WORK WITH RAW PROXIMITY ANALYSIS OUTPUT!
"""

class Oligos():

	""" Class containing properties of oligonucleotides. 
		The oligonucleotide information comes from Tongjia's proximity output
	"""

	def __init__(self,
			     oligo_name,
			     oligo_orientation,
			     start_location,
			     end_location,
			     subject_id,
			     subject_name,
			     subject_accession,
			     query_percent,
			     e_value):

		self.oligo_name = oligo_name
		self.oligo_orientation = oligo_orientation
		self.start_location = start_location
		self.end_location = end_location
		self.subject_id = subject_id
		self.subject_name = subject_name
		self.subject_accession = subject_accession
		self.query_percent = query_percent
		self.e_value = e_value

class Matches():

	""" Main purpose is to match a unique subject id to the oligonucleotide pair
		(reverse) and (forward) that it corresponds to. The subject id comes from
		an automated generated number by GenomeQuest (www.gqlifescience.com).
		The oligonucleotide information comes from Tongjia's proximity output.

		Class containing the matches that are flagged as matched when
		the string 'RESULT MATCHED' is found anywhere in the data file.
	"""

	def __init__(self, subject_id, oligo_pairs):
		self.subject_id = subject_id
		self.oligo_pairs = oligo_pairs

class ProximityRip():

	""" Main Program that allows execution and writing 
		of the output into a more friendly readable context

		Main Methods:
			(1) EXECUTE
			(2) WRITE_TO_FILE
	"""

	def __init__(self):
		self.match_list = {}

	def execute(self, filename):

		oligo_list = self.__create_oligos(filename)
		self.match_list = self.__create_matches(filename, oligo_list)

	def write_to_file(self, save_as):

		with open(save_as, "w") as save_file:

			headers = ["Subject ID", "Accession", "Description", 
			           "Oligo 1 (FOR)", "Oligo 2(REV)", "Oligo 1(E-val)",
			           "Oligo 1(Query Percent)", "Oligo 1(Start Location)",
			           "Oligo 1(End Location", "Oligo 2(E-val)", "Oligo 2(Query Percent)",
			           "Oligo 2(Start Location", "Oligo 2(End Location"]

			save_file.write("\t".join(headers))
			save_file.write("\n")

			for subjects, values in self.match_list.iteritems():
				output = [subjects]
				# General Properties | Could've used either probe pairs
				output.append(values.oligo_pairs[0].subject_accession)
				output.append(values.oligo_pairs[0].subject_name)
				output.append(values.oligo_pairs[0].oligo_name)
				output.append(values.oligo_pairs[1].oligo_name)
				# Oligo 1 Properties
				output.append(values.oligo_pairs[0].e_value)
				output.append(values.oligo_pairs[0].query_percent)
				output.append(values.oligo_pairs[0].start_location)
				output.append(values.oligo_pairs[0].end_location)
				# Oligo 2 Properties
				output.append(values.oligo_pairs[1].e_value)
				output.append(values.oligo_pairs[1].query_percent)
				output.append(values.oligo_pairs[1].start_location)
				output.append(values.oligo_pairs[1].end_location)

				final_output = "\t".join(output)
				save_file.write(final_output)
				save_file.write("\n")

		print "Proximity Rip has been completed."
		print "File has been saved as %s" % save_as			

	def __create_oligos(self, filename):

		"""
		Creates oligonucleotide objects with properties from the
		proximity analysis output file. File must be in a .tsv format
		with the following data format:

		General File Format of an output. {x} notate
		that there is data available, but we choose to ignore it.

		Line # | Data Type
		0) SUBJECT ID
		1) ACCESSION
		2) E-VALUE
 		3) {x}
		4) QUERY ID PERCENT
		5) {x}
		6) {x}
		7) {x}
		8) {x}
		9) OLIGO NAME
		10) SUBJECT NAME
		11) START LOCATION
		12) END LOCATION

		"""

		import csv

		oligo_obj_list = {}

		with open(filename) as q:

			csv_reader = csv.reader(q, delimiter='\t')

			for lines in csv_reader:

				if lines:
					if not (lines[1].find('too far') > -1 or 
						    lines[1].find('RESULT') > -1 or 
						    lines[1].find('reversed') > -1):
			

						oligo_orientation = lines[9].split("_")[0]

						oligo_factory = Oligos(lines[9],					# OLIGO NAME
											   oligo_orientation,			# OLIGO ORIENTATION
											   lines[11],					# SUBJECT START LOCATION
											   lines[12],					# SUBJECT END LOCATION
											   lines[0],					# SUBJECT ID
											   lines[10],					# SUBJECT NAME
											   lines[1],					# SUBJECT ACCESSION
											   lines[4],					# QUERY ID
											   lines[2])					# E VALUE

						if lines[0] in oligo_obj_list:
							oligo_obj_list[lines[0]].append(oligo_factory)
						else:
							oligo_obj_list[lines[0]] = [oligo_factory]

		return oligo_obj_list 

	def __create_matches(self, filename, list_of_oligos):
		""" Gets a list of subject ID for MATCHES
		were found. 'Matches' are determined by the keywords
		'RESULT MATCHED'.

		Args:
			filename - proximity analysis file <str>
			list_of_olgios - dict containing subj id (key) & list of oligos (value) <dict>
		Returns:
			a dictionary with subject id (key) and corresponding oligo matches (value) <dict>

		"""

		import csv

		list_of_result_matches = {}

		with open(filename) as m:
			csv_reader = csv.reader(m, delimiter='\t')
			for lines in csv_reader:
				if lines:
					if (lines[1].find('RESULT MATCHED') > -1):
						curr_subject_id = lines[0]
						start_locations = self.__parse_locations(lines[1])
						oligo_pairs = self.__oligo_pairs(start_locations, list_of_oligos[curr_subject_id])
						list_of_result_matches[curr_subject_id] = Matches(curr_subject_id, oligo_pairs)

		return list_of_result_matches

	def __parse_locations(self, result_string):
		""" Used for the purpose of matching back to which
			oligo object the result is referring to 

			ex: FOR_OLIGOA:3400 REV_OLIGOB:3900
			returns (3400,3900)

			Args:
				result_string - the string that contains 'RESULT MATCHED'
			Returns:
				tuple of locations for forward and reverse oligonucleotides (for,rev) <tuple>
		"""
		split_result = result_string.split(" ")

		forward_oligo_name = split_result[4]
		reverse_oligo_name = split_result[6]

		for_location_index = forward_oligo_name.index(":")

		rev_location_index = reverse_oligo_name.index(":")
		rev_end_location_index = reverse_oligo_name.index("*")

		forward_start_location = forward_oligo_name[for_location_index+1:]
		reverse_start_location = reverse_oligo_name[rev_location_index+1:rev_end_location_index]

		return [forward_start_location, reverse_start_location]

	def __oligo_pairs(self, locations, oligo_list):
		""" Get the correct oligos from a list of locations 

		Args:
			oligo_list - oligo list for a specific subject id <list>
			locations - oligo start locations that were previously parsed <list>
		Returns:
			only the oligos that have a matching start location <tuple>
		"""

		FORWARD = ""
		REVERSE = ""

		for oligo in oligo_list:
			if oligo.start_location in locations:
				if oligo.oligo_orientation == "FOR":
					FORWARD = oligo
				else:
					REVERSE = oligo

		pair_of_oligos = [FORWARD,REVERSE]

		return pair_of_oligos

if __name__ == '__main__':

	from sys import argv

	print ""
	print "########## PROXIMITY RIP 1.0 ##########"
	print "Internal Use Only"
	print "Requirements: Filtered output of proximity analysis"
	print ""
	print "FLOW OF PROCESS:"
	print "ProximityAnalysis.jar --> ProximityFilter2.py --> ProximityRip.py"
	print "#######################################"
	print "Use: python proximityrip.py <filtered proximity file> <save as>"
	print ""

	p = ProximityRip()

	try:
		p.execute(argv[1])
		p.write_to_file(argv[2])
	except IOError:
		print "File does not exist."
	except IndexError:
		print "ERROR: Incorrect use of program. Check parameters."
		print "Use: python proximityrip.py <filtered proximity file> <save as>"