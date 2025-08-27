class fragmenter:
    """Class for fragmenting molecules based on predefined SMARTS patterns and algorithms.

    This class provides methods to perform fragmentation on molecules using different algorithms
    ('simple', 'complete', or 'combined'). It utilizes RDKit for chemical structure handling and supports
    deep copying of objects, adjacency matrix calculations, and recursive searches for complete fragmentations.

    Attributes:
        algorithm (str): The fragmentation algorithm to use.
        n_max_fragmentations_to_find (int): Maximum number of fragmentations to find.
        n_heavy_atoms_cuttoff (int): Heavy atom cutoff for fragmentation.
        match_hydrogens (bool): Whether to consider hydrogens during matching.
        fragmentation_scheme (dict): Mapping of group numbers to SMARTS patterns or lists of SMARTS patterns.
        function_to_choose_fragmentation (callable): Function used to select the best fragmentation from multiple options.
        fragmentation_scheme_order (list): Order in which to search for fragmentation groups.
        _fragmentation_scheme_group_number_lookup (dict): Lookup mapping SMARTS patterns to group numbers.
        _fragmentation_scheme_pattern_lookup (dict): Lookup mapping SMARTS patterns to RDKit Mol objects.
        _adjacency_matrix_cache (dict): Cache for storing computed adjacency matrices keyed by canonical SMILES.
    """

    from rdkit import Chem
    import marshal as marshal
    from rdkit.Chem import rdmolops
    import warnings

    @staticmethod
    def deep_copy(obj):
        """Perform a deep copy of an object using marshal serialization.

        This method is based on the approach described in:
        https://stackoverflow.com/questions/24756712/deepcopy-is-extremely-slow

        Args:
            obj: The object to deep copy.

        Returns:
            A deep copy of the input object.
        """
        return fragmenter.marshal.loads(fragmenter.marshal.dumps(obj))

    @staticmethod
    def get_heavy_atom_count(mol):
        """Count the number of heavy (non-hydrogen) atoms in a molecule.

        Args:
            mol: An RDKit molecule object.

        Returns:
            int: The count of heavy atoms in the molecule.
        """
        return sum(atom.GetAtomicNum() != 1 for atom in mol.GetAtoms())

    def get_substruct_matches(
        self,
        mol_searched_for,
        mol_searched_in,
        canonical_SMILES_searched_in,
        atom_indices_to_which_new_matches_have_to_be_adjacent,
    ):
        """Retrieve substructure matches from a molecule that are adjacent to specified atom indices.

        This method performs a substructure search and then filters the matches to include only those that are
        adjacent to the provided set of atom indices.

        Args:
            mol_searched_for: An RDKit Mol representing the substructure pattern to search for.
            mol_searched_in: An RDKit Mol in which to search for the substructure.
            canonical_SMILES_searched_in (str): The canonical SMILES representation of the molecule (used for caching).
            atom_indices_to_which_new_matches_have_to_be_adjacent (set): Set of atom indices; new matches must be adjacent
                to at least one of these atoms.

        Returns:
            list: A list of tuples, where each tuple contains the atom indices corresponding to a matching substructure.
                  Returns an empty list if no matches are found or if the molecule is too small.
        """
        if mol_searched_in.GetNumAtoms() < mol_searched_for.GetNumAtoms():
            return []

        matches = mol_searched_in.GetSubstructMatches(mol_searched_for)
        if not matches:
            return []

        if not atom_indices_to_which_new_matches_have_to_be_adjacent:
            return list(matches)

        valid_matches = []
        adjancency_matrix = self.get_adjacency_matrix(
            mol_searched_in, canonical_SMILES_searched_in
        )
        for match in matches:
            for i in match:
                if adjancency_matrix[i].intersection(
                    atom_indices_to_which_new_matches_have_to_be_adjacent
                ):
                    valid_matches.append(match)
                    break

        return valid_matches

    def __init__(
        self,
        fragmentation_scheme={},
        fragmentation_scheme_order=None,
        match_hydrogens=False,
        algorithm="",
        n_heavy_atoms_cuttoff=-1,
        function_to_choose_fragmentation=False,
        n_max_fragmentations_to_find=-1,
    ):
        """Initialize the fragmenter with a fragmentation scheme and algorithm parameters.

        Args:
            fragmentation_scheme (dict): Dictionary mapping group numbers (int) to SMARTS strings or list of SMARTS strings.
            fragmentation_scheme_order (list, optional): Order in which to search for fragmentation groups.
                Defaults to None.
            match_hydrogens (bool, optional): Whether to include hydrogens in matching. Defaults to False.
            algorithm (str): Algorithm to use for fragmentation; must be one of 'simple', 'complete', or 'combined'.
            n_heavy_atoms_cuttoff (int, optional): Maximum number of heavy atoms allowed for fragmentation in
                'complete' or 'combined' algorithms. Defaults to -1.
            function_to_choose_fragmentation (callable or bool, optional): Function to select the best fragmentation
                among possible fragmentations. Required for 'combined' algorithms. Defaults to False.
            n_max_fragmentations_to_find (int, optional): Maximum number of fragmentations to find.
                Defaults to -1 (no limit) for the 'simple' algorithm.

        Raises:
            TypeError: If fragmentation_scheme is not a dictionary.
            ValueError: If fragmentation_scheme is empty.
            ValueError: If algorithm is not one of 'simple', 'complete', or 'combined'.
            ValueError: If n_max_fragmentations_to_find is set for algorithm 'simple'.
            ValueError: If n_heavy_atoms_cuttoff or function_to_choose_fragmentation is not properly set for
                'complete' or 'combined' algorithms.
            TypeError: If function_to_choose_fragmentation is not callable or returns an incorrect type.
        """
        if not type(fragmentation_scheme) is dict:
            raise TypeError(
                "fragmentation_scheme must be a dctionary with integers as keys and either strings or list of strings as values."
            )

        if len(fragmentation_scheme) == 0:
            raise ValueError("fragmentation_scheme must be provided.")

        if not algorithm in ["simple", "complete", "combined"]:
            raise ValueError("Algorithm must be either simple ,complete or combined.")

        if algorithm == "simple":
            if n_max_fragmentations_to_find != -1:
                raise ValueError(
                    "Setting n_max_fragmentations_to_find only makes sense with complete or combined algorithm."
                )

        self.algorithm = algorithm

        if algorithm in ["combined", "complete"]:
            if n_heavy_atoms_cuttoff == -1:
                raise ValueError(
                    "n_atoms_cuttoff needs to be specified for complete or combined algorithms."
                )
            if algorithm == "combined":
                if function_to_choose_fragmentation == False:
                    raise ValueError(
                        "function_to_choose_fragmentation needs to be specified for complete or combined algorithms."
                    )

                if not callable(function_to_choose_fragmentation):
                    raise TypeError(
                        "function_to_choose_fragmentation needs to be a function."
                    )
                else:
                    if type(function_to_choose_fragmentation([{}, {}])) not in [
                        dict,
                        list,
                    ]:
                        raise TypeError(
                            "function_to_choose_fragmentation needs to take a list of fragmentations and return one fragmentation or a list of fragmentations."
                        )

            if n_max_fragmentations_to_find != -1:
                if n_max_fragmentations_to_find < 1:
                    raise ValueError(
                        "n_max_fragmentations_to_find has to be 1 or higher."
                    )

        if fragmentation_scheme_order is None:

            # create automagic fragmentation_scheme_order from sizes of groups from largest to smallest and specificity of SMARTS
            scheme_descriptors = []
            for i, SMARTS in fragmentation_scheme.items():
                if isinstance(SMARTS, list):
                    SMARTS = SMARTS[0]
                mol_SMARTS = fragmenter.Chem.MolFromSmarts(SMARTS)
                n_atoms_SMARTS = 0
                for atom in mol_SMARTS.GetAtoms():
                    try:
                        n_atoms_SMARTS += 1
                        atom_query = atom.DescribeQuery()
                        if "AtomHCount" not in atom_query:
                            continue
                        n_Hs = int(atom_query.split("AtomHCount")[1].split("=")[0])
                        n_atoms_SMARTS += n_Hs
                    except Exception:
                        pass
                scheme_descriptors.append((i, n_atoms_SMARTS, len(SMARTS)))

            fragmentation_scheme_order = [
                i
                for i, p1, p2 in sorted(
                    scheme_descriptors, key=lambda x: (x[1], x[2]), reverse=True
                )
            ]

            if algorithm in ["simple", "combined"]:
                self.warnings.warn(
                    "No especific fragmentation_scheme_order was given, groups were sorted by group size from largest to smallest, you might get better results if you specify the order in which the groups are searched for."
                )

        self.n_max_fragmentations_to_find = n_max_fragmentations_to_find

        self.n_heavy_atoms_cuttoff = n_heavy_atoms_cuttoff

        self.match_hydrogens = match_hydrogens

        self.fragmentation_scheme = fragmentation_scheme

        self.function_to_choose_fragmentation = function_to_choose_fragmentation

        # create lookup dictionaries to faster finding a group number
        self._fragmentation_scheme_group_number_lookup = {}
        self._fragmentation_scheme_pattern_lookup = {}
        self.fragmentation_scheme_order = fragmentation_scheme_order
        self._adjacency_matrix_cache = {}

        for group_number, list_SMARTS in fragmentation_scheme.items():

            if type(list_SMARTS) is not list:
                list_SMARTS = [list_SMARTS]

            for SMARTS in list_SMARTS:
                if SMARTS != "":
                    self._fragmentation_scheme_group_number_lookup[SMARTS] = (
                        group_number
                    )

                    mol_SMARTS = fragmenter.Chem.MolFromSmarts(SMARTS)
                    self._fragmentation_scheme_pattern_lookup[SMARTS] = mol_SMARTS

    def fragment(self, SMILES_or_molecule):
        """Fragment a molecule using the configured fragmentation scheme and algorithm.

        This method processes the molecule—converting from a SMILES string if necessary—and applies the
        fragmentation algorithm to each disconnected fragment.

        Args:
            SMILES_or_molecule (str or Mol): The molecule to fragment, provided as a SMILES string or an RDKit Mol object.

        Returns:
            tuple: A tuple containing:
                - fragmentation (dict): Mapping of group numbers to counts of matches found.
                - success (bool): True if the fragmentation was successful (i.e., all atoms were assigned), False otherwise.
                - fragmentation_matches (dict): Mapping of group numbers to lists of match tuples (atom index sequences).

        Raises:
            ValueError: If a provided SMILES string is invalid.
        """
        if isinstance(SMILES_or_molecule, str):
            complete_mol = fragmenter.Chem.MolFromSmiles(SMILES_or_molecule)
            complete_mol = (
                fragmenter.Chem.AddHs(complete_mol)
                if self.match_hydrogens
                else complete_mol
            )
            is_valid_SMILES = complete_mol is not None

            if not is_valid_SMILES:
                raise ValueError("Following SMILES is not valid: " + SMILES_or_molecule)

        else:
            complete_mol = SMILES_or_molecule

        # iterate over all separated molecules
        success = []
        fragmentation = {}
        fragmentation_matches = {}
        frags = list(fragmenter.rdmolops.GetMolFrags(complete_mol, asMols=True))
        for mol in frags:
            SMILES = fragmenter.Chem.MolToSmiles(mol)
            this_mol_fragmentation, this_mol_success = self.__get_fragmentation(
                mol, SMILES
            )

            for SMARTS, matches in this_mol_fragmentation.items():
                group_number = self._fragmentation_scheme_group_number_lookup[SMARTS]

                if not group_number in fragmentation:
                    fragmentation[group_number] = 0
                    fragmentation_matches[group_number] = []

                fragmentation[group_number] += len(matches)
                fragmentation_matches[group_number].extend(matches)

            success.append(this_mol_success)

        return fragmentation, all(success), fragmentation_matches

    def fragment_complete(self, SMILES_or_molecule):
        """Perform complete fragmentation on a single-fragment molecule.

        This method only accepts molecules with a single connected component and attempts to find
        all possible complete fragmentations that cover the entire molecule.

        Args:
            SMILES_or_molecule (str or Mol): The molecule to fragment, provided as a SMILES string or an RDKit Mol object.

        Returns:
            tuple: A tuple containing:
                - fragmentations (list): A list of dictionaries mapping group numbers to counts of matches in each complete fragmentation.
                - success (bool): True if at least one complete fragmentation was found, False otherwise.
                - fragmentations_matches (list): A list of dictionaries mapping group numbers to lists of match tuples for each fragmentation.

        Raises:
            ValueError: If a provided SMILES string is invalid or if the molecule consists of multiple fragments.
        """
        if type(SMILES_or_molecule) is str:
            mol_SMILES = fragmenter.Chem.MolFromSmiles(SMILES_or_molecule)
            mol_SMILES = (
                fragmenter.Chem.AddHs(mol_SMILES)
                if self.match_hydrogens
                else mol_SMILES
            )
            is_valid_SMILES = mol_SMILES is not None

            if not is_valid_SMILES:
                raise ValueError("Following SMILES is not valid: " + SMILES_or_molecule)

        else:
            mol_SMILES = SMILES_or_molecule

        if len(fragmenter.rdmolops.GetMolFrags(mol_SMILES)) != 1:
            raise ValueError(
                "fragment_complete does not accept multifragment molecules."
            )

        temp_fragmentations, success = self.__complete_fragmentation(
            mol_SMILES, fragmenter.Chem.MolToSmiles(mol_SMILES)
        )

        fragmentations = []
        fragmentations_matches = []
        for temp_fragmentation in temp_fragmentations:
            fragmentation = {}
            fragmentation_matches = {}
            for SMARTS, matches in temp_fragmentation.items():
                group_number = self._fragmentation_scheme_group_number_lookup[SMARTS]

                fragmentation[group_number] = len(matches)
                fragmentation_matches[group_number] = matches

            fragmentations.append(fragmentation)
            fragmentations_matches.append(fragmentation_matches)

        return fragmentations, success, fragmentations_matches

    def __get_fragmentation(self, mol, canonical_SMILES):
        """Retrieve a fragmentation solution for a molecule.

        This method first attempts a simple fragmentation and, if unsuccessful, resorts to a complete fragmentation.
        If multiple fragmentations are found, the provided function selects the preferred one.

        Args:
            mol: An RDKit Mol object representing the molecule to fragment.
            canonical_SMILES (str): The canonical SMILES string of the molecule.

        Returns:
            tuple: A tuple containing:
                - fragmentation (dict): Mapping of SMARTS patterns to lists of match tuples.
                - success (bool): True if a complete fragmentation was achieved, False otherwise.
        """
        success = False
        fragmentation = {}
        if self.algorithm in ["simple", "combined"]:
            fragmentation, success = self.__simple_fragmentation(mol, canonical_SMILES)

        if success:
            return fragmentation, success

        if self.algorithm in ["combined", "complete"]:
            fragmentations, success = self.__complete_fragmentation(
                mol, canonical_SMILES
            )

            if success and self.algorithm == "combined":
                fragmentation = self.function_to_choose_fragmentation(fragmentations)

        return fragmentation, success

    def __simple_fragmentation(self, mol, canonical_SMILES):
        """Attempt to fragment a molecule using the simple algorithm.

        The simple algorithm searches for non-overlapping substructure matches iteratively and,
        if necessary, cleans the molecule to improve coverage.

        Args:
            mol: An RDKit Mol object representing the molecule.
            canonical_SMILES (str): The canonical SMILES string of the molecule.

        Returns:
            tuple: A tuple containing:
                - fragmentation (dict): Mapping of SMARTS patterns to lists of match tuples.
                - success (bool): True if the entire molecule was successfully fragmented, False otherwise.
        """
        if self.match_hydrogens:
            target_atom_count = mol.GetNumAtoms()
        else:
            target_atom_count = fragmenter.get_heavy_atom_count(mol)

        success = False
        fragmentation = {}

        fragmentation, atom_indices_included_in_fragmentation = (
            self.__search_non_overlapping_solution(
                mol, canonical_SMILES, {}, set(), set(), {}
            )
        )
        success = len(atom_indices_included_in_fragmentation) == target_atom_count

        # if not successful, clean up molecule and search again
        level = 1
        while not success:
            cleaned_fragmentation, cleaned_atom_indices = (
                self.__clean_molecule_surrounding_unmatched_atoms(
                    mol,
                    canonical_SMILES,
                    fragmentation,
                    atom_indices_included_in_fragmentation,
                    level,
                )
            )

            if not cleaned_atom_indices:
                break

            cleaned_fragmentation, cleaned_atom_indices = (
                self.__search_non_overlapping_solution(
                    mol,
                    canonical_SMILES,
                    cleaned_fragmentation,
                    cleaned_atom_indices,
                    cleaned_atom_indices,
                    {},
                )
            )

            if len(cleaned_atom_indices) == target_atom_count:
                success = True
                fragmentation = cleaned_fragmentation
            else:
                level += 1

        if not success and fragmentation:
            # in case of incomplete fragmentation
            # try one match at a time whether a specific match
            # leads to an incomplete fragmentation
            groups_leading_to_incomplete_fragmentation = {}
            for SMARTS in fragmentation:
                for match in fragmentation[SMARTS]:
                    groups_leading_to_incomplete_fragmentation.setdefault(
                        SMARTS, []
                    ).append(match)

                    cleaned_fragmentation, cleaned_atom_indices = (
                        self.__search_non_overlapping_solution(
                            mol,
                            canonical_SMILES,
                            cleaned_fragmentation,
                            cleaned_atom_indices,
                            cleaned_atom_indices,
                            groups_leading_to_incomplete_fragmentation,
                        )
                    )

                    if len(cleaned_atom_indices) == target_atom_count:
                        success = True
                        fragmentation = cleaned_fragmentation
                        break

                if success:
                    break

        return fragmentation, success

    def __search_non_overlapping_solution(
        self,
        mol_searched_in,
        canonical_SMILES_searched_in,
        fragmentation,
        atom_indices_included_in_fragmentation,
        atom_indices_to_which_new_matches_have_to_be_adjacent,
        groups_leading_to_incomplete_fragmentations,
    ):
        """Iteratively search for non-overlapping substructure matches to fragment the molecule.

        This method updates the provided fragmentation mapping and set of assigned atom indices by
        iteratively searching for additional matches that do not overlap with previously assigned atoms.

        Args:
            mol_searched_in: The RDKit Mol object in which to search.
            canonical_SMILES_searched_in (str): The canonical SMILES representation of the molecule.
            fragmentation (dict): Current mapping of SMARTS patterns to match lists.
            atom_indices_included_in_fragmentation (set): Set of atom indices already included in the fragmentation.
            atom_indices_to_which_new_matches_have_to_be_adjacent (set): Set of atom indices that new matches must be adjacent to.
            groups_leading_to_incomplete_fragmentations (dict): keys are SMARTS, vlues are a list of tuples of found atom indices.
        Returns:
            tuple: A tuple containing:
                - fragmentation (dict): Updated fragmentation mapping.
                - atom_indices_included_in_fragmentation (set): Updated set of atom indices included in the fragmentation.
        """
        prev_count = -1

        while prev_count != len(atom_indices_included_in_fragmentation):
            prev_count = len(atom_indices_included_in_fragmentation)

            for group_number in self.fragmentation_scheme_order:
                list_SMARTS = self.fragmentation_scheme[group_number]

                if isinstance(list_SMARTS, str):
                    list_SMARTS = [list_SMARTS]

                for SMARTS in list_SMARTS:
                    if SMARTS:
                        fragmentation, atom_indices_included_in_fragmentation = (
                            self.__get_next_non_overlapping_match(
                                mol_searched_in,
                                canonical_SMILES_searched_in,
                                SMARTS,
                                fragmentation,
                                atom_indices_included_in_fragmentation,
                                atom_indices_to_which_new_matches_have_to_be_adjacent,
                                groups_leading_to_incomplete_fragmentations,
                            )
                        )

        return fragmentation, atom_indices_included_in_fragmentation

    def __get_next_non_overlapping_match(
        self,
        mol_searched_in,
        canonical_SMILES_searched_in,
        SMARTS,
        fragmentation,
        atom_indices_included_in_fragmentation,
        atom_indices_to_which_new_matches_have_to_be_adjacent,
        groups_leading_to_incomplete_fragmentations,
    ):
        """Find and record the next non-overlapping substructure match for a given SMARTS pattern.

        Searches for substructure matches in the molecule that are adjacent to specified atom indices and
        do not overlap with atoms already assigned in the fragmentation. If a valid match is found, it is added
        to the fragmentation mapping.

        Args:
            mol_searched_in: The RDKit Mol object to search within.
            canonical_SMILES_searched_in (str): The canonical SMILES representation of the molecule.
            SMARTS (str): The SMARTS pattern to search for.
            fragmentation (dict): Current mapping of SMARTS patterns to lists of matches.
            atom_indices_included_in_fragmentation (set): Set of atom indices already assigned.
            atom_indices_to_which_new_matches_have_to_be_adjacent (set): Set of atom indices that new matches must be adjacent to.
            groups_leading_to_incomplete_fragmentations (dict): keys are SMARTS, vlues are a list of tuples of found atom indices.

        Returns:
            tuple: A tuple containing:
                - fragmentation (dict): Updated fragmentation mapping with the new match, if found.
                - atom_indices_included_in_fragmentation (set): Updated set of atom indices assigned.
        """
        mol_searched_for = self._fragmentation_scheme_pattern_lookup[SMARTS]

        matches = self.get_substruct_matches(
            mol_searched_for,
            mol_searched_in,
            canonical_SMILES_searched_in,
            atom_indices_to_which_new_matches_have_to_be_adjacent or None,
        )

        if not matches:
            return fragmentation, atom_indices_included_in_fragmentation

        for match in matches:
            if atom_indices_included_in_fragmentation.isdisjoint(match):
                if (
                    SMARTS not in groups_leading_to_incomplete_fragmentations
                    or match not in groups_leading_to_incomplete_fragmentations[SMARTS]
                ):
                    fragmentation.setdefault(SMARTS, []).append(match)
                    atom_indices_included_in_fragmentation.update(match)

        return fragmentation, atom_indices_included_in_fragmentation

    def get_adjacency_matrix(self, mol, canonical_SMILES):
        """Compute and cache the adjacency matrix for a molecule.

        The adjacency matrix is represented as a list of sets, where each set contains the indices of neighboring atoms.

        Args:
            mol: An RDKit Mol object for which to compute the adjacency matrix.
            canonical_SMILES (str): The canonical SMILES representation of the molecule used as a cache key.

        Returns:
            list: A list of sets, each representing the neighboring atom indices for each atom in the molecule.
        """
        if canonical_SMILES in self._adjacency_matrix_cache:
            return self._adjacency_matrix_cache[canonical_SMILES]

        adjancency_matrix = []
        for atom in mol.GetAtoms():
            adjancency_matrix.append(
                {neighbor.GetIdx() for neighbor in atom.GetNeighbors()}
            )

        self._adjacency_matrix_cache[canonical_SMILES] = adjancency_matrix
        return adjancency_matrix

    def __clean_molecule_surrounding_unmatched_atoms(
        self,
        mol_searched_in,
        canonical_SMILES_searched_in,
        fragmentation,
        atom_indices_included_in_fragmentation,
        level,
    ):
        """Clean the fragmentation by removing matches that are adjacent to unmatched atoms.

        This method refines the current fragmentation by iteratively removing substructure matches
        that border atoms which have not been assigned any match. The cleaning is performed for the specified
        number of levels.

        Args:
            mol_searched_in: The RDKit Mol object being fragmented.
            canonical_SMILES_searched_in (str): The canonical SMILES representation of the molecule.
            fragmentation (dict): Current mapping of SMARTS patterns to match lists.
            atom_indices_included_in_fragmentation (set): Set of atom indices already included in the fragmentation.
            level (int): The number of cleaning iterations to perform.

        Returns:
            tuple: A tuple containing:
                - fragmentation (dict): The cleaned fragmentation mapping.
                - atom_indices_included_in_fragmentation (set): Updated set of atom indices included after cleaning.
        """
        total_atoms = fragmenter.get_heavy_atom_count(mol_searched_in)

        for _ in range(level):

            atoms_missing = (
                set(range(total_atoms)) - atom_indices_included_in_fragmentation
            )

            new_fragmentation = fragmenter.deep_copy(fragmentation)

            atom_to_smart_mapping = {}
            for smart, atoms_found in fragmentation.items():
                for atoms in atoms_found:
                    for atom in atoms:
                        atom_to_smart_mapping.setdefault(atom, []).append(
                            (smart, atoms)
                        )

            adjancency_matrix = self.get_adjacency_matrix(
                mol_searched_in, canonical_SMILES_searched_in
            )
            for atom_index in atoms_missing:
                for neighbor_index in adjancency_matrix[atom_index]:
                    if neighbor_index in atom_to_smart_mapping:
                        for smart, atoms in atom_to_smart_mapping[neighbor_index]:
                            if (
                                smart in new_fragmentation
                                and atoms in new_fragmentation[smart]
                            ):
                                new_fragmentation[smart].remove(atoms)

            new_fragmentation = {k: v for k, v in new_fragmentation.items() if v}

            atom_indices_included_in_fragmentation = {
                atom
                for atoms_list in new_fragmentation.values()
                for atoms in atoms_list
                for atom in atoms
            }

            fragmentation = new_fragmentation

        return fragmentation, atom_indices_included_in_fragmentation

    def __complete_fragmentation(self, mol, canonical_SMILES):
        """Perform a complete fragmentation search on a molecule within a heavy atom count cutoff.

        This method recursively searches for all possible complete fragmentations that cover the entire molecule,
        subject to a heavy atom cutoff.

        Args:
            mol: An RDKit Mol object representing the molecule.
            canonical_SMILES (str): The canonical SMILES representation of the molecule.

        Returns:
            tuple: A tuple containing:
                - completed_fragmentations (list): A list of fragmentation mappings (SMARTS to match lists).
                - success (bool): True if at least one complete fragmentation was found, False otherwise.
        """
        heavy_atom_count = fragmenter.get_heavy_atom_count(mol)

        if heavy_atom_count > self.n_heavy_atoms_cuttoff:
            return {}, False

        if self.match_hydrogens:
            target_atom_count = mol.GetNumAtoms()
        else:
            target_atom_count = heavy_atom_count

        completed_fragmentations = []
        groups_leading_to_incomplete_fragmentations = []
        (
            completed_fragmentations,
            groups_leading_to_incomplete_fragmentations,
            incomplete_fragmentation_found,
        ) = self.__get_next_non_overlapping_adjacent_match_recursively(
            mol,
            canonical_SMILES,
            target_atom_count,
            completed_fragmentations,
            groups_leading_to_incomplete_fragmentations,
            {},
            set(),
            set(),
            self.n_max_fragmentations_to_find,
        )
        success = len(completed_fragmentations) > 0

        return completed_fragmentations, success

    def __get_next_non_overlapping_adjacent_match_recursively(
        self,
        mol_searched_in,
        canonical_SMILES_searched_in,
        target_atom_count,
        completed_fragmentations,
        groups_leading_to_incomplete_fragmentations,
        fragmentation_so_far,
        atom_indices_included_in_fragmentation_so_far,
        atom_indices_to_which_new_matches_have_to_be_adjacent,
        n_max_fragmentations_to_find=-1,
    ):
        """Recursively search for non-overlapping adjacent substructure matches to achieve complete fragmentation.

        This method explores possible fragmentations by recursively adding valid substructure matches until the entire
        molecule is covered. It also tracks incomplete fragmentations to avoid redundant searches.

        Args:
            mol_searched_in: The RDKit Mol object to fragment.
            canonical_SMILES_searched_in (str): The canonical SMILES representation of the molecule.
            target_atom_count (int): Total number of atoms (or heavy atoms) that must be covered.
            completed_fragmentations (list): List of complete fragmentation mappings found so far.
            groups_leading_to_incomplete_fragmentations (list): List of fragmentation groups that led to incomplete fragmentations.
            fragmentation_so_far (dict): Current partial fragmentation mapping.
            atom_indices_included_in_fragmentation_so_far (set): Set of atom indices already included in the partial fragmentation.
            atom_indices_to_which_new_matches_have_to_be_adjacent (set): Set of atom indices that new matches must be adjacent to.
            n_max_fragmentations_to_find (int, optional): Maximum number of complete fragmentations to find.
                Defaults to -1 (no limit).

        Returns:
            tuple: A tuple containing:
                - completed_fragmentations (list): Updated list of complete fragmentation mappings.
                - groups_leading_to_incomplete_fragmentations (list): Updated list of groups leading to incomplete fragmentations.
                - incomplete_fragmentation_found (bool): True if an incomplete fragmentation was detected, False otherwise.
        """
        n_completed_fragmentations = len(completed_fragmentations)
        incomplete_fragmentation_found = False
        complete_fragmentation_found = False

        if len(completed_fragmentations) == n_max_fragmentations_to_find:
            return (
                completed_fragmentations,
                groups_leading_to_incomplete_fragmentations,
                incomplete_fragmentation_found,
            )

        for group_number in self.fragmentation_scheme_order:
            list_SMARTS = self.fragmentation_scheme[group_number]

            if complete_fragmentation_found:
                break

            if type(list_SMARTS) is not list:
                list_SMARTS = [list_SMARTS]

            for SMARTS in list_SMARTS:
                if complete_fragmentation_found:
                    break

                if SMARTS != "":
                    matches = self.get_substruct_matches(
                        self._fragmentation_scheme_pattern_lookup[SMARTS],
                        mol_searched_in,
                        canonical_SMILES_searched_in,
                        atom_indices_included_in_fragmentation_so_far,
                    )

                    for match in matches:

                        # only allow non-overlapping matches
                        # as get_substruct_matches does not test this
                        all_atoms_are_unassigned = (
                            atom_indices_included_in_fragmentation_so_far.isdisjoint(
                                match
                            )
                        )
                        if not all_atoms_are_unassigned:
                            continue

                        # only allow matches that do not contain groups leading to incomplete matches
                        for (
                            groups_leading_to_incomplete_fragmentation
                        ) in groups_leading_to_incomplete_fragmentations:
                            if self.__is_fragmentation_subset_of_other_fragmentation(
                                groups_leading_to_incomplete_fragmentation,
                                fragmentation_so_far,
                            ):
                                return (
                                    completed_fragmentations,
                                    groups_leading_to_incomplete_fragmentations,
                                    incomplete_fragmentation_found,
                                )

                        # only allow matches that will lead to new fragmentations
                        use_this_match = True
                        n_found_groups_so_far = len(fragmentation_so_far)

                        for completed_fragmentation in completed_fragmentations:

                            if not SMARTS in completed_fragmentation:
                                continue

                            if n_found_groups_so_far == 0:
                                use_this_match = (
                                    not self.__is_match_contained_in_fragmentation(
                                        match, SMARTS, completed_fragmentation
                                    )
                                )
                            else:
                                if self.__is_fragmentation_subset_of_other_fragmentation(
                                    fragmentation_so_far, completed_fragmentation
                                ):
                                    use_this_match = (
                                        not self.__is_match_contained_in_fragmentation(
                                            match, SMARTS, completed_fragmentation
                                        )
                                    )

                            if not use_this_match:
                                break

                        if not use_this_match:
                            continue

                        # make a deepcopy here, otherwise the variables are modified down the road
                        this_SMARTS_fragmentation_so_far = fragmenter.deep_copy(
                            fragmentation_so_far
                        )
                        this_SMARTS_atom_indices_included_in_fragmentation_so_far = (
                            atom_indices_included_in_fragmentation_so_far.copy()
                        )

                        if not SMARTS in this_SMARTS_fragmentation_so_far:
                            this_SMARTS_fragmentation_so_far[SMARTS] = []

                        this_SMARTS_fragmentation_so_far[SMARTS].append(match)
                        this_SMARTS_atom_indices_included_in_fragmentation_so_far.update(
                            match
                        )

                        # only allow matches that do not contain groups leading to incomplete matches
                        for (
                            groups_leading_to_incomplete_match
                        ) in groups_leading_to_incomplete_fragmentations:
                            if self.__is_fragmentation_subset_of_other_fragmentation(
                                groups_leading_to_incomplete_match,
                                this_SMARTS_fragmentation_so_far,
                            ):
                                use_this_match = False
                                break

                        if not use_this_match:
                            continue

                        # if the complete molecule has not been fragmented, continue to do so
                        if (
                            len(
                                this_SMARTS_atom_indices_included_in_fragmentation_so_far
                            )
                            < target_atom_count
                        ):
                            (
                                completed_fragmentations,
                                groups_leading_to_incomplete_fragmentations,
                                incomplete_fragmentation_found,
                            ) = self.__get_next_non_overlapping_adjacent_match_recursively(
                                mol_searched_in,
                                canonical_SMILES_searched_in,
                                target_atom_count,
                                completed_fragmentations,
                                groups_leading_to_incomplete_fragmentations,
                                this_SMARTS_fragmentation_so_far,
                                this_SMARTS_atom_indices_included_in_fragmentation_so_far,
                                this_SMARTS_atom_indices_included_in_fragmentation_so_far,
                                n_max_fragmentations_to_find,
                            )
                            break

                        # if the complete molecule has been fragmented, save and return
                        if (
                            len(
                                this_SMARTS_atom_indices_included_in_fragmentation_so_far
                            )
                            == target_atom_count
                        ):
                            completed_fragmentations.append(
                                this_SMARTS_fragmentation_so_far
                            )
                            complete_fragmentation_found = True
                            break

        # if until here no new fragmentation was found check whether an incomplete fragmentation was found
        if n_completed_fragmentations == len(completed_fragmentations):

            if not incomplete_fragmentation_found:

                incomplete_matched_groups = {}

                if len(atom_indices_included_in_fragmentation_so_far) > 0:
                    unassignes_atom_indices = set(
                        range(0, target_atom_count)
                    ).difference(atom_indices_included_in_fragmentation_so_far)
                    adjancency_matrix = self.get_adjacency_matrix(
                        mol_searched_in, canonical_SMILES_searched_in
                    )
                    for atom_index in unassignes_atom_indices:
                        for neighbor_atom_index in adjancency_matrix[atom_index]:
                            for (
                                found_smarts,
                                found_matches,
                            ) in fragmentation_so_far.items():
                                for found_match in found_matches:
                                    if neighbor_atom_index in found_match:
                                        if (
                                            not found_smarts
                                            in incomplete_matched_groups
                                        ):
                                            incomplete_matched_groups[found_smarts] = []

                                        if (
                                            found_match
                                            not in incomplete_matched_groups[
                                                found_smarts
                                            ]
                                        ):
                                            incomplete_matched_groups[
                                                found_smarts
                                            ].append(found_match)

                    is_subset_of_groups_already_found = False
                    indexes_to_remove = []

                    for group_index, groups_leading_to_incomplete_match in enumerate(
                        groups_leading_to_incomplete_fragmentations
                    ):
                        is_subset_of_groups_already_found = (
                            self.__is_fragmentation_subset_of_other_fragmentation(
                                incomplete_matched_groups,
                                groups_leading_to_incomplete_match,
                            )
                        )
                        if is_subset_of_groups_already_found:
                            indexes_to_remove.append(group_index)

                    for index in sorted(indexes_to_remove, reverse=True):
                        del groups_leading_to_incomplete_fragmentations[index]

                    groups_leading_to_incomplete_fragmentations.append(
                        incomplete_matched_groups
                    )
                    groups_leading_to_incomplete_fragmentations = sorted(
                        groups_leading_to_incomplete_fragmentations, key=len
                    )

                    incomplete_fragmentation_found = True

        return (
            completed_fragmentations,
            groups_leading_to_incomplete_fragmentations,
            incomplete_fragmentation_found,
        )

    def __is_fragmentation_subset_of_other_fragmentation(
        self, fragmentation, other_fragmentation
    ):
        """Determine if one fragmentation mapping is a subset of another.

        This method checks whether every match in the given fragmentation exists within the other fragmentation.

        Args:
            fragmentation (dict): A fragmentation mapping to test.
            other_fragmentation (dict): Another fragmentation mapping to compare against.

        Returns:
            bool: True if 'fragmentation' is a subset of 'other_fragmentation', False otherwise.
        """
        if not fragmentation:
            return False

        if len(other_fragmentation) < len(fragmentation):
            return False

        for found_SMARTS, matches in fragmentation.items():
            if found_SMARTS not in other_fragmentation:
                return False

            found_matches_set = {frozenset(i) for i in matches}
            found_other_matches_set = {
                frozenset(i) for i in other_fragmentation[found_SMARTS]
            }

            if not found_matches_set.issubset(found_other_matches_set):
                return False

        return True

    def __is_match_contained_in_fragmentation(self, match, SMARTS, fragmentation):
        """Check if a specific match is already contained within a fragmentation mapping for a given SMARTS pattern.

        Args:
            match (tuple): A tuple of atom indices representing the substructure match.
            SMARTS (str): The SMARTS pattern associated with the match.
            fragmentation (dict): The fragmentation mapping containing SMARTS keys and their corresponding matches.

        Returns:
            bool: True if the match is already present in the fragmentation mapping, False otherwise.
        """
        if SMARTS not in fragmentation:
            return False

        found_matches_set = {frozenset(i) for i in fragmentation[SMARTS]}

        return frozenset(match) in found_matches_set


if __name__ == "__main__":
    main()
