using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace VinaCarbResults
{
    class Program
    {
		public const int linesAddedByResults = 3;

		static bool saveBestPDBQT = true;
		static bool saveBestModes = true;
		static bool saveAllModes = true;

		static string enz = "";			//d

        static void Main(string[] args)
        {
            #region Args population
			// Paths: 3klkG.enzyme	3klk.enzyme	3hz3.enzyme
            var allEnzymeInfo = Directory.GetFiles(Directory.GetCurrentDirectory(), "*.enzyme", SearchOption.TopDirectoryOnly)
				.Select(x => new FileInfo(x));

			// Paths: /ENumbers	/CatLib
			var possibleDirsInfo = Directory.GetDirectories(Directory.GetCurrentDirectory(), "*", SearchOption.TopDirectoryOnly)
				.Except(allEnzymeInfo.Select(x => x.FullName.Substring(0, x.FullName.Length - x.Extension.Length)))
				.Select(x => new DirectoryInfo(x));

			var dirs = PromptArg(args, 0, "What substrate?", possibleDirsInfo.Select(x => x.Name).ToArray(), "");

			if (dirs != "all")
				possibleDirsInfo = possibleDirsInfo.Where(x => x.Name == dirs);

			// Paths: /Enumbers/E100Curcumin	/CatLib/CG
			var possiblePosesInfo = new List<DirectoryInfo>();
			foreach (var dir in possibleDirsInfo)
			{
				var possiblePoses = Directory.GetDirectories(dir.FullName, "*", SearchOption.TopDirectoryOnly)
					.Select(x => new DirectoryInfo(x));

				possiblePosesInfo.AddRange(possiblePoses);
			}

			if (!possiblePosesInfo.Any())
				return;	//TODO: Error out?

			var poses = PromptArg(args, 2, "What pose?", possiblePosesInfo.Select(x => x.Name).ToArray(), "");

			// Remove any poses the user didn't want
			if (poses != "all")
				possiblePosesInfo.RemoveAll(x => x.Name != poses);

			// Loop through remaining poses to see what enzymes were used to dock them
			var knownResults = new List<string>();
			foreach (var pose in possiblePosesInfo)
			{
				knownResults.AddRange(
					Directory.GetFiles(pose.FullName, "*.pdbqt", SearchOption.TopDirectoryOnly)
						.Select(x => new FileInfo(x).Name)
					);
			}

			enz = PromptArg(args, 1, "What enzyme?", allEnzymeInfo
				.Where(x => knownResults.Any(y => NameEnz(y).enzymeName + ".enzyme" == x.Name))
				.Select(x => Path.GetFileNameWithoutExtension(x.Name))
				.ToArray(), "");

			//Remove all poses that were not run for any enzymes we're interested in
			if (enz != "all")
			{
				possiblePosesInfo.RemoveAll(x =>
				{
					return !Directory.GetFiles(x.FullName, "*.pdbqt", SearchOption.TopDirectoryOnly)
						.Any(y => NameEnz(new FileInfo(y).Name).enzymeName == enz);
				});

				enz += ".enzyme";
			}
			#endregion

			/* x x x x x x   ITERATE ALL PDBQT AND LOG   x x x x x x x */
			//Run all poses one-by-one (print pose summary of each)
			// Paths: /Enumbers/E100Curcumin	/CatLib/CG
			possiblePosesInfo.ForEach(x => ProcessPose(x));

			//ExtractAngle(enzyme);	Not implemented

			/* x x x x x x   WRITE DATA TO CONSOLE AND TO FILE   x x x x x x x */
			possiblePosesInfo.ForEach(x => PerformSaves(x));

			//Stall exit until Esc or Backspace if not started with args
			if (!args.Any())
				StallTillEscBsp();
		}

		static void PerformSaves(DirectoryInfo pose)
		{
			var poseInst = GetPose(pose.FullName);
			if (saveBestPDBQT)
			{
				foreach (var kvp in poseInst.EnzymeModes)
				{
					var builder = new StringBuilder();
					int i = 1;

					//Two best binding modes CANNOT be the same run and bindingNr!
					//This is because only one atom can be in the active site - so the same [run, bindNr] is never two best modes
					//So no duplicates!

					foreach (var mode in poseInst.BestModes(kvp.Key))
                    {
						var resultFileName = kvp.Key + (mode.run > 1 ? "_" + mode.run : "") + ".pdbqt";

						//Open each file and extract lines relevant to mode.run
						//lines from (mode.bindingNr - 1) * linePerResult to (mode.bindingNr * linePerResult)

						//Add "MODE _" line in order
						builder.AppendLine("MODE " + i);
						i++;

						//Read lines of bindingNr (MODE _)
						File.ReadLines(pose.FullName + "//" + resultFileName)
							.Skip((mode.bindingNr - 1) * poseInst.linePerResult + 1)	//Skip all previous modes and ignore MODE _ line
							.Take(poseInst.linePerResult - 1)	//Take all lines except for the MODE line
							.ToList().ForEach(x => builder.AppendLine(x));	//Append all lines to StringBuilder
					}

					//TODO: Pick a nicer file (format) for results, and make sure this code ignores such files
					File.WriteAllText(pose.FullName + "//_" + kvp.Key + ".pdb", builder.ToString());
				}
			}

			if (saveBestModes)
			{
				var builder = new StringBuilder();
				poseInst.BestModes(ref builder);
				File.WriteAllText(pose.FullName + "//_best.tsv", builder.ToString());
			}

			if (saveAllModes)
			{
				var builder = new StringBuilder();
				poseInst.AllModes(ref builder);
				File.WriteAllText(pose.FullName + "//_all.tsv", builder.ToString());
			}
        }

		static void ProcessPose(DirectoryInfo pose)
		{
            //Console.WriteLine();
            Console.Write(pose.Name);

            #region Populate results
            var possibleResults = Directory.GetFiles(pose.FullName,
					(enz == "all") ? "*.pdbqt" : enz+"*.pdbqt",		//First filter: enz*.pdbqt
					SearchOption.TopDirectoryOnly)
				.Select(x => new FileInfo(x));

			//Make sure the file starts with the correct enzyme name
			if (enz != "all")
				possibleResults = possibleResults.Where(x => enz == NameEnz(x.Name).enzymeName);
			//x => x.StartsWith(enz + ".") || x.StartsWith(enz + "_");
			#endregion

			var poseInst = GetPose(pose.FullName);
			var atomsList = poseInst.atoms.Select(x => x.lineNrPDBQT + 2).ToList(); //Two lines are added in results, so shift indices
			var totalAtoms = atomsList.Count;
			var entryLength = poseInst.linePerResult;

			foreach (var result in possibleResults)
			{
				var split = NameEnz(result.Name);

				var enzymeSphere = activeSite(split.enzymeName);
				if (enzymeSphere.sqRadius == 0)
					return;   //The enzyme file didn't exist or doesn't have a radius comparer -- no reason to check result

				//Copy ALL log file lines
				var simpleModes = File.ReadLines(Regex.Replace(result.FullName, ".pdbqt", ".log"))
					.SkipWhile(x => !x.StartsWith("--"))
					.Skip(1)                                    //Skip to the part where values are listed
					.TakeWhile(x => !x.StartsWith("Writing"))   //Take until the end
					.Select(x =>
					{
						var part = x.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
						return new Mode(split.run, int.Parse(part[0]), float.Parse(part[3]), float.Parse(part[1]));
					}).ToList();

				if (!simpleModes.Any())
					continue;

				File.ReadLines(result.FullName)
					.Where((val, index) => atomsList.Contains(index % entryLength))     //Only read if the atom is known
																						// NOTE!! The index now has a different meaning!
					.Select((val, index) => new { dist = Vector3.SqDistance(LinePos(val), enzymeSphere.position), index })
					.Where(x => x.dist <= enzymeSphere.sqRadius)
					.All(x =>
						{
							var mode = x.index / totalAtoms; //+ (x.index % entryLength == 0 ? 0 : 1);
						var atomNr = x.index % totalAtoms;  //Continuous integer index of atomsList
						simpleModes[mode] = new Mode(simpleModes[mode], atomNr);
							return true;
						});

				//Finally add all modes to the pose for this enzyme
				if (!poseInst.EnzymeModes.ContainsKey(split.enzymeName))
					poseInst.EnzymeModes.Add(split.enzymeName, simpleModes);
				else
					poseInst.EnzymeModes[split.enzymeName].AddRange(simpleModes);

				//Console.Write("|");
			}

			//Console.WriteLine();
			Console.WriteLine(GetPose(pose.FullName).ToString());
		}

		//Try on: http://regexstorm.net/tester !!
		static Regex regex1 = new Regex("(.+?)(?:_([0-9]+))?(?:\\.pdbqt)");

		static (string enzymeName, int run) NameEnz(string name)
		{
			var match = regex1.Match(name);
			if (!match.Success)
				throw new ArgumentException(name + " is not of the right form");
			var run = (match.Groups[2].Value != "") ? int.Parse(match.Groups[2].Value) : 1;
			return (match.Groups[1].Value, run);
		}

		static Vector3 LinePos(string line)
		{
			var part = line.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
			if (int.TryParse(part[3], out _))
			{
				//	(HETERO)ATOM	ATOMNR	NAME	GROUPNR		X	Y	Z
				//	0				1		2		3			4	5	6
				return new Vector3(float.Parse(part[4]), float.Parse(part[5]), float.Parse(part[6]));
			}
			else if (int.TryParse(part[4], out _))
			{
				//	(HETERO)ATOM	ATOMNR	NAME	NAME	GROUPNR		X	Y	Z
				//	0				1		2		3		4			5	6	7
				return new Vector3(float.Parse(part[5]), float.Parse(part[6]), float.Parse(part[7]));
			}

			//	(HETERO)ATOM	ATOMNR	NAME	NAME	AA	GROUPNR		X	Y	Z
			//	0				1		2		3		4	5			6	7	8
			return new Vector3(float.Parse(part[6]), float.Parse(part[7]), float.Parse(part[8]));
		}

		static AtomProps LineProps(string line, int index)
		{
			var part = line.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

			var atomNr = int.Parse(part[1]);
			var atomName = part[2];

			if (int.TryParse(part[3], out int groupNr2))
			{
				//	(HETERO)ATOM	ATOMNR	NAME	GROUPNR
				//	0				1		2		3
				return new AtomProps(index, atomNr, atomName, groupNr2, "UNK");
			}
			else if (int.TryParse(part[4], out int groupNr))
			{
				//	(HETERO)ATOM	ATOMNR	NAME	NAME	GROUPNR
				//	0				1		2		3		4
				return new AtomProps(index, atomNr, atomName, groupNr, part[3]);
			}

			//	(HETERO)ATOM	ATOMNR	NAME	NAME	NAME	GROUPNR
			//	0				1		2		3		4		5
			return new AtomProps(index, atomNr, atomName, int.Parse(part[5]), part[3]);
		}

		static void StallTillEscBsp()
		{
			var key = Console.ReadKey();
			while (key.Key != ConsoleKey.Backspace && key.Key != ConsoleKey.Escape)
			{
				key = Console.ReadKey();
			}
		}

		//static Dictionary<string, List<AtomProps>> poseAtoms = new Dictionary<string, List<AtomProps>>(); //Key: dir; Value: all atoms
		static Dictionary<string, Sphere> activeSiteSphere = new Dictionary<string, Sphere>();
		static Dictionary<string, Pose> dirPoses = new Dictionary<string, Pose>(); //Key: DirectoryInfo.FullName; Value: Pose

		static Pose GetPose(string path)
		{
			if (!dirPoses.ContainsKey(path))
			{
				var p = new Pose();

				var newPath = path + ".pdbqt";
				bool simplePDBQT = File.Exists(newPath);

				//Try to extract atom information from a different .pdbqt (e.g from the results), by counting distance between models and extracting oxygens.
				if (!simplePDBQT)
					newPath = Directory.GetFiles(path, "*.pdbqt", SearchOption.TopDirectoryOnly).FirstOrDefault();

				//.pdbqt doesn't exist (anymore)! NOR DO RESULTS
				if (newPath == default)
				{
					//TODO: Store error
					return null;
				}

				var lines = File.ReadLines(newPath);

				if (!simplePDBQT)
					lines = lines.Skip(2).TakeWhile(x => !x.StartsWith("ENDMDL"));

				p.lineNumberPerEntryPDBQT = lines.Count();

				p.atoms.AddRange(lines       //Read filepath
					.Select((value, index) => new { index, value })     //Store line number before "Where"
					.Where(x => x.value.EndsWith("O") || x.value.EndsWith("OA"))    //Select only oxygens
					.Select(x => LineProps(x.value, x.index)));

				dirPoses.Add(path, p);
			}
			return dirPoses[path];
		}

		/*
		//Order list to correct order for header
		//Olist = Olist.OrderBy(x => x.groupNr+"."+x.atomNr).ToList();

		dirPoses.Add(filePath, p);
		poseAtoms.Add(filePath, p.atoms);
		*/

		static Sphere activeSite(string enzyme)
        {
			if (!activeSiteSphere.ContainsKey(enzyme))
			{
				if (File.Exists(enzyme + ".enzyme"))
				{

					var radiusLine = File.ReadLines(enzyme + ".enzyme").First(x => x.Contains("--radius"));
					if (radiusLine == null)
						activeSiteSphere.Add(enzyme, default);

					var splits = radiusLine.Split(' ');
					var radius = float.Parse(splits[1]);
					var pos = new Vector3(float.Parse(splits[2]), float.Parse(splits[3]), float.Parse(splits[4]));
					activeSiteSphere.Add(enzyme, new Sphere(pos, radius));
				}
				else
				{
					//TODO: Store error


					//Set to default sphere
					activeSiteSphere.Add(enzyme, default);
				}
			}
			return activeSiteSphere[enzyme];
        }

		static string PromptArg(string[] args, int index, string prompt, string[] options, string defaultVal, bool addAll = true)
		{
			var testVal = defaultVal;
			if (args.Length > index)
            {
				testVal = args[index].ToLower();
				if (testVal == "all")
					return testVal;
            }

			if (testVal.Length == 0 || (!options.Contains(testVal) && testVal != "all"))
			{
				if (addAll)
					options = (new[] { "all" }).Concat(options).ToArray();

				testVal = AutoPrompt.PromptForInput(prompt + " [all or autocomplete] :", options);
			}
			return testVal;
		}

    }
}
