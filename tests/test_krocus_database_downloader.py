import unittest
import os
import tempfile
import logging
from unittest.mock import Mock, patch, MagicMock
from krocus.KrocusDatabaseDownloader import KrocusDatabaseDownloader

class TestKrocusDatabaseDownloader(unittest.TestCase):

    def test_initialization(self):
        """Test KrocusDatabaseDownloader object initialization"""
        options = Mock()
        options.list_species = False
        options.species = None
        options.output_directory = 'test_output'
        options.verbose = False
        
        downloader = KrocusDatabaseDownloader(options)
        self.assertEqual(downloader.output_directory, 'test_output')
        self.assertFalse(downloader.list_species)

    def test_verbose_logging(self):
        """Test that verbose mode sets DEBUG logging level"""
        options = Mock()
        options.list_species = False
        options.species = None
        options.output_directory = 'test_output'
        options.verbose = True
        
        downloader = KrocusDatabaseDownloader(options)
        self.assertEqual(downloader.logger.level, logging.DEBUG)

    def test_non_verbose_logging(self):
        """Test that non-verbose mode sets ERROR logging level"""
        options = Mock()
        options.list_species = False
        options.species = None
        options.output_directory = 'test_output'
        options.verbose = False
        
        downloader = KrocusDatabaseDownloader(options)
        self.assertEqual(downloader.logger.level, logging.ERROR)

    @patch('krocus.KrocusDatabaseDownloader.PubmlstGetter')
    def test_run_list_species(self, mock_pubmlst):
        """Test run method with list_species option"""
        options = Mock()
        options.list_species = True
        options.species = None
        options.output_directory = 'test_output'
        options.verbose = False
        
        mock_getter = MagicMock()
        mock_pubmlst.return_value = mock_getter
        
        downloader = KrocusDatabaseDownloader(options)
        downloader.run()
        
        mock_getter.print_available_species.assert_called_once()

    @patch('krocus.KrocusDatabaseDownloader.PubmlstGetter')
    def test_run_download_species(self, mock_pubmlst):
        """Test run method with species download"""
        options = Mock()
        options.list_species = False
        options.species = 'Salmonella enterica'
        options.output_directory = 'test_output'
        options.verbose = False
        
        mock_getter = MagicMock()
        mock_pubmlst.return_value = mock_getter
        
        downloader = KrocusDatabaseDownloader(options)
        downloader.run()
        
        mock_getter.get_species_files.assert_called_once_with('Salmonella enterica', 'test_output')

    @patch('krocus.KrocusDatabaseDownloader.PubmlstGetter')
    def test_run_no_valid_options(self, mock_pubmlst):
        """Test run method with no valid options"""
        options = Mock()
        options.list_species = False
        options.species = None
        options.output_directory = 'test_output'
        options.verbose = False
        
        mock_getter = MagicMock()
        mock_pubmlst.return_value = mock_getter
        
        downloader = KrocusDatabaseDownloader(options)
        downloader.run()
        
        # Should not call any methods on pubmlst getter
        mock_getter.print_available_species.assert_not_called()
        mock_getter.get_species_files.assert_not_called()

if __name__ == '__main__':
    unittest.main()
