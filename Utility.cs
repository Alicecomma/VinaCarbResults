using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace VinaCarbResults
{

	public struct Vector3
	{
		float x;
		float y;
		float z;

		public Vector3(float x, float y, float z)
		{
			this.x = x;
			this.y = y;
			this.z = z;
		}

		public static float SqDistance(Vector3 a, Vector3 b)
		{
			float diff_x = a.x - b.x;
			float diff_y = a.y - b.y;
			float diff_z = a.z - b.z;
			return (float)(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z);
		}

		public static float Distance(Vector3 a, Vector3 b)
		{
			float diff_x = a.x - b.x;
			float diff_y = a.y - b.y;
			float diff_z = a.z - b.z;
			return (float)Math.Sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z);
		}

		public float sqrMagnitude => x * x + y * y + z * z;
	}

	public struct Sphere
	{
		public Vector3 position;
		float radius;
		public float sqRadius;

		public Sphere(Vector3 pos, float rad)
		{
			position = pos;
			radius = rad;
			sqRadius = (float)Math.Pow(rad, 2);
		}
	}

	public class Pose
	{
		public const int linesAddedByResults = 3;

		public int lineNumberPerEntryPDBQT;
		public int linePerResult => lineNumberPerEntryPDBQT + linesAddedByResults;

		//public Dictionary<string, List<AtomModes>> AtomsForEnzyme;  //Key: enzyme name; Val: list of atoms

		public List<AtomProps> atoms = new List<AtomProps>();   //Links [i] to atomNr, atomName, groupNr, groupName

		/// <summary>
		/// Key: 3klk	3hz3	3klkG;	Value: Mode
		/// </summary>
		public Dictionary<string, List<Mode>> EnzymeModes = new Dictionary<string, List<Mode>>();

		public List<Mode> BestModes(string enzyme)
		{
			var list = new List<Mode>(atoms.Count);
			if (EnzymeModes.TryGetValue(enzyme, out var modes))
			{
				for (int i = 0; i < atoms.Count; i++)
				{
					var knownAtoms = modes.Where(x => x.atomIndex == i); //Ascending order so lowest is first

					if (knownAtoms.Count() > 0)
						list.Add(knownAtoms.OrderBy(x => x.carbScore).First());
					else
						list.Add(default);
				}
			}
			return list;
        }

		public override string ToString()
        {
			var t = "\t";
			var builder = new StringBuilder();

			BestModes(ref builder);
			//AllModes(ref builder);

			return builder.ToString();
		}

		public void AllModes(ref StringBuilder builder)
		{
			var t = "\t";

			if (!EnzymeModes.Any())
			{
				builder.AppendLine("No results");
				return;
			}

			/* x x x x x x   CREATE BANNER   x x x x x x x */

			builder.AppendLine("Enzyme" + t + "Label" + t + "Name" + t + "Run" + t + "Mode" + t + "VC" + t + "AV");
				/*	Enzyme	Label	Name	Run	Mode	VC		AV
				 *			(nr)	(id)	nr	nr		kcal/mol */
			builder.AppendLine(t + "(nr)" + t + "(id)" + t + "nr" + t + "nr" + t + "kcal/mol");

			/* x x x x x x   LIST OF ALL MODES   x x x x x x x */
			foreach (var kvp in EnzymeModes)
			{
				builder.AppendLine().Append(kvp.Key);

				foreach (var mode in kvp.Value.Where(x => x.atomIndex > -1).OrderBy(x => x.carbScore))
				{
					var atom = atoms[mode.atomIndex];
					builder.AppendLine(t + atom.ToString() + t + mode.ToString());

					/*	3klkG	GLC1(1)	O1(16)	1	16	-7.5	-8.2
					 *			GLC2(2)	O6(13)	1	15	-6.0	-8.0
					 *			.(.)	.(.)	.	.	.	.
					 */
				}
			}
		}

		public void BestModes(ref StringBuilder builder)
		{
			var t = "\t";

			/* x x x x x x   CREATE THE BANNER   x x x x x x x */
			//List<string> columns = new List<string>(atoms.Count);
			//foreach (var atom in atoms.GroupBy(x => x.groupNr+"("+x.atomNr)).ToList())
			//{

			if (atoms.Count < 1)
			{
				builder.AppendLine("No oxygens");
				return;
			}

			string line1 = string.Empty;
			string line2 = string.Empty;

			for (int i = 0; i < atoms.Count; i++)
			{
				/*
				 *		group(label)					|	group(label)
				 *		id(name)	id(name)	id(name)	|	id(name)	id(name)
				 */

				line1 += t + atoms[i].groupName + "(" + atoms[i].groupNr + ")";
				line2 += t + atoms[i].atomName + "(" + atoms[i].atomNr + ")";

				// TODO: Make optimal ordering of atom indices for display (group labels)
			}

			builder.AppendLine(line1).AppendLine(line2);

			/* x x x x x x   BEST BINDING MODES   x x x x x x x */
			foreach (var kvp in EnzymeModes)
			{
				string line3 = kvp.Key;
				string line4 = string.Empty;
				string line5 = string.Empty;
				string line6 = string.Empty;

				/*
				 *	Enz	-7.0	-6.5	-3.0	-4.5	kcal/mol
				 *		1	3	16	10	nr
				 *		1	2	1	3	run
				 */

				bool[] p = new bool[atoms.Count];   //Store whether VinaCarb score != AutoDock Vina score

				var bests = BestModes(kvp.Key);
				for (int i = 0; i < atoms.Count; i++)
				{
					if (bests[i].bindingNr != default) //&& bests[i].run != default)
					{
						line3 += t + bests[i].carbScore;
						p[i] = bests[i].carbScore != bests[i].vinaScore;
						line4 += t + (string)(p[i] ? bests[i].vinaScore.ToString() : "");
						line5 += t + bests[i].bindingNr;
						line6 += t + bests[i].run;
					}
					else
					{
						line3 += t;
						line4 += t;
						line5 += t;
						line6 += t;
					}
				}

				var anyVC = p.Any(x => x);

				builder.AppendLine().AppendLine(line3 + t + "kcal/mol" + (anyVC ? " (VC)" : ""));
				if (anyVC) builder.AppendLine(line4 + t + "kcal/mol (AV)");
				builder.AppendLine(line5 + t + "mode").AppendLine(line6 + t + "run");

			}
		}
		
	}

	public struct AtomProps
	{
		public int lineNrPDBQT;
		public int atomNr;
		public string atomName;
		public int groupNr;
		public string groupName;

		public AtomProps(int lineNrPDBQT, int atomNr, string atomName, int groupNr, string groupName)
        {
			this.lineNrPDBQT = lineNrPDBQT;
			this.atomNr = atomNr;
			this.atomName = atomName;
			this.groupNr = groupNr;
			this.groupName = groupName;
        }

        public override string ToString()
        {
			return groupName + "(" + groupNr + ")\t" + atomName + "("+ atomNr + ")";
        }
    }

	public class AtomModes
    {
		public AtomProps props;
		public List<Mode> modes;

		public Mode BestMode => modes.OrderBy(x => x.carbScore).LastOrDefault();
    }

	public struct Mode
	{
		public int run;         //If this value is zero, then the Mode has never been instantiated
		public int bindingNr;	//If this value is zero, then the Mode has never been instantiated
		public float vinaScore;
		public float carbScore;
		public int atomIndex;

		public Mode (int run, int bindingNr, float vinaScore, float carbScore, int atomIndex = -1)
        {
			this.run = run;
			this.bindingNr = bindingNr;
			this.vinaScore = vinaScore;
			this.carbScore = carbScore;
			this.atomIndex = atomIndex;
        }

		public Mode (Mode mode, int atomIndex)
        {
			this.run = mode.run;
			this.bindingNr = mode.bindingNr;
			this.vinaScore = mode.vinaScore;
			this.carbScore = mode.carbScore;
			this.atomIndex = atomIndex;
        }

        public override string ToString()
        {
			return run + "\t" + bindingNr + "\t" + carbScore + (carbScore != vinaScore ? "\t" + vinaScore : "");
        }
    }
}
